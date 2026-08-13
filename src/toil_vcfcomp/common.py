"""
Shared helpers for running the comparison tools (hap.py, truvari, rtg vcfeval, aardvark) either
in docker, in singularity, or as plain binaries on the PATH.

Commands are always written with container paths: the work directory is bind-mounted at
CONTAINER_MNT, so a command refers to its inputs as <CONTAINER_MNT>/<basename>.  In local mode
run_tool() rewrites those back to real paths, so the callers in vcfcomp.py do not have to know
which mode is in effect.

This module deliberately does not import toil at the top level, so that the standalone `vcfcomp`
entry point works in an environment that has no toil installed.
"""

import os, sys
import subprocess
import errno
import hashlib
import shutil
import tempfile
import time
from datetime import datetime
from urllib.parse import urlparse

try:
    from toil.realtimeLogger import RealtimeLogger
except ImportError:
    RealtimeLogger = None

# Where the work directory is bind-mounted inside the container.  /mnt rather than /data because
# only the aardvark image actually ships a /data directory, and a singularity or apptainer install
# without overlay/underlay support silently SKIPS a bind whose destination does not exist -- the
# tool then reports its inputs as missing, for a reason nothing in the log explains.  /mnt exists
# in all four tool images, which is why cactus binds there too.
CONTAINER_MNT = '/mnt'
MNT = CONTAINER_MNT + '/'


def log_info(msg, max_len=1500):
    """ log to toil's realtime logger when we are inside a job, otherwise to stderr.  truncate
    long messages so they are less likely to get dropped (same trick as cactus_realtime_log) """
    if len(msg) > max_len:
        msg = msg[:max_len - 207] + ' <...> ' + msg[-200:]
    stamped = '{}: {}'.format(datetime.now(), msg)
    if RealtimeLogger:
        RealtimeLogger.info(stamped)
    else:
        sys.stderr.write(stamped + '\n')


def makeURL(path_or_url):
    """ toil wants a scheme on everything it imports or exports """
    if urlparse(path_or_url).scheme == '':
        return 'file://' + os.path.abspath(path_or_url)
    else:
        return path_or_url


# --- binaries mode ---

def add_binaries_options(parser):
    """ the one option block shared by `vcfcomp` and `toil-vcfcomp` """
    parser.add_argument('--binariesMode', choices=['docker', 'local', 'singularity'], default=None,
                        help='how to run hap.py/truvari/rtg/aardvark. autodetected if not given: '
                        'docker if the docker daemon answers, else singularity, else local')


def setup_binaries(options):
    """ resolve --binariesMode once, in main(), and stash it in the environment so the job
    functions can get at it without carrying options around.

    Autodetection runs wherever main() runs, which under a cluster batch system is the submit
    node -- not where the jobs will actually execute.  A submit node commonly has docker (it is
    where images get built) while the compute nodes have only singularity, so preferring docker
    there would pick a runtime no worker can use.  Prefer singularity whenever the work is going
    somewhere else.
    """
    mode = getattr(options, 'binariesMode', None)
    if mode is None:
        mode = os.environ.get('VCFCOMP_BINARIES_MODE', None)
    if mode is None:
        batch_system = str(getattr(options, 'batchSystem', 'single_machine')).lower()
        distributed = batch_system not in ['single_machine', 'singlemachine']
        has_singularity = bool(shutil.which('singularity') or shutil.which('apptainer'))
        if distributed and has_singularity:
            mode = 'singularity'
        elif _tool_works(['docker', 'info']):
            mode = 'docker'
        elif has_singularity:
            mode = 'singularity'
        else:
            mode = 'local'
        log_info('Autodetected --binariesMode {}'.format(mode))
    # keep it on options too, so a --restart cannot silently resolve to a different runtime than
    # the run it is resuming
    if hasattr(options, 'binariesMode'):
        options.binariesMode = mode
    os.environ['VCFCOMP_BINARIES_MODE'] = mode
    return mode


def binaries_mode():
    return os.environ.get('VCFCOMP_BINARIES_MODE', 'docker')


def _tool_works(cmd):
    # timeout because `docker info` hangs rather than failing when the daemon socket is present
    # but unresponsive, which would stall startup indefinitely
    try:
        subprocess.check_call(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, timeout=30)
        return True
    except (subprocess.CalledProcessError, subprocess.TimeoutExpired, FileNotFoundError, OSError):
        return False


def _singularity_exe():
    return shutil.which('singularity') or shutil.which('apptainer') or 'singularity'


# --- running tools ---

def run_tool(parameters, image=None, work_dir=None, check_output=False, outfile=None,
             job_memory=None, file_store=None):
    """ run a tool command.  `parameters` is an argv list whose paths are written as
    <CONTAINER_MNT>/<name>; work_dir is what gets mounted there.  `image` is the docker specifier,
    ignored in local mode.  returns stdout when check_output is set. """
    mode = binaries_mode()
    if work_dir is None:
        work_dir = os.getcwd()
    work_dir = os.path.abspath(work_dir)

    if mode == 'docker':
        if not image:
            raise RuntimeError('--binariesMode docker needs a docker image for: {}'.format(parameters[0]))
        # note: --interactive but never -t.  a tty cannot be allocated under a batch executor,
        # and asking for one fails with "the input device is not a TTY"
        call = ['docker', 'run', '--rm', '--interactive', '--log-driver=none',
                '-u', '{}:{}'.format(os.getuid(), os.getgid()),
                '-v', '{}:{}'.format(work_dir, CONTAINER_MNT), '-w', CONTAINER_MNT,
                image] + parameters
    elif mode == 'singularity':
        if not image:
            raise RuntimeError('--binariesMode singularity needs a docker image for: {}'.format(parameters[0]))
        sandbox = singularity_sandbox(image, work_dir=work_dir, file_store=file_store)
        call = [_singularity_exe(), '--silent', 'exec', '-u',
                '-B', '{}:{}'.format(work_dir, CONTAINER_MNT), '--pwd', CONTAINER_MNT,
                sandbox] + parameters
    else:
        assert mode == 'local'
        call = [_localize_path(p, work_dir) for p in parameters]

    log_info('Running: "{}"'.format(' '.join(call)))
    start_time = time.time()

    if check_output:
        result = subprocess.check_output(call, encoding='utf-8')
    elif outfile:
        # stage under a temporary name and rename only once the command has exited cleanly, so
        # a truncated output can never masquerade as a finished one
        outfile_tmp = outfile + '.part'
        with open(outfile_tmp, 'w') as out_fh:
            subprocess.check_call(call, stdout=out_fh)
        os.rename(outfile_tmp, outfile)
        result = None
    else:
        subprocess.check_call(call)
        result = None

    run_time = time.time() - start_time
    msg = 'Successfully ran: "{}" in {} seconds'.format(' '.join(call), round(run_time, 4))
    if job_memory:
        msg += ' with job-memory {}'.format(job_memory)
    log_info(msg)
    return result


def _localize_path(param, work_dir):
    """ /data/x -> <work_dir>/x, for running outside a container """
    if param.startswith(MNT):
        return os.path.join(work_dir, param[len(MNT):])
    return param


def singularity_sandbox(image, work_dir=None, file_store=None, cache_dir=None):
    """ build (or find) a singularity sandbox directory for a docker image.

    Concurrent `singularity pull`s that share a cache directory corrupt each other
    (https://github.com/sylabs/singularity/issues/3634), which is exactly what happens when a few
    hundred slurm jobs start at once.  So we keep our own cache, build into a temp directory with
    a private SINGULARITY_CACHEDIR, and atomically rename it into place.  Sandbox directories
    rather than .sif files, since only those work with user namespaces.  Cribbed from cactus,
    which cribbed it from toil-vg. """
    if cache_dir is None:
        cache_dir = os.path.join(os.environ.get('SINGULARITY_CACHEDIR',
                                                os.path.join(os.path.expanduser('~'), '.singularity')),
                                 'toil-vcfcomp')
    os.makedirs(cache_dir, exist_ok=True)

    spec = image
    if not spec.startswith('/') and '://' not in spec:
        spec = 'docker://' + spec

    sandbox_dir = os.path.join(cache_dir, '{}.sandbox'.format(hashlib.sha256(spec.encode('utf-8')).hexdigest()))
    if os.path.exists(sandbox_dir):
        return sandbox_dir

    # 'building.' prefix so a leftover is identifiable as debris rather than a real cache entry
    temp_sandbox_dir = tempfile.mkdtemp(dir=cache_dir, prefix='building.')
    download_env = os.environ.copy()
    if file_store:
        download_env['SINGULARITY_CACHEDIR'] = file_store.getLocalTempDir()
    else:
        download_env['SINGULARITY_CACHEDIR'] = tempfile.mkdtemp(dir=work_dir if work_dir else cache_dir)

    build_cmd = [_singularity_exe(), 'build', '-s', '-F', temp_sandbox_dir, spec]
    log_info('Building singularity sandbox: "{}"'.format(' '.join(build_cmd)))
    start_time = time.time()
    try:
        subprocess.check_call(build_cmd, env=download_env)
    except BaseException:
        # a build killed by a time limit or a registry hiccup otherwise leaves a multi-GB partial
        # sandbox in a shared $HOME cache under a name nothing will ever clean up
        shutil.rmtree(temp_sandbox_dir, ignore_errors=True)
        shutil.rmtree(download_env['SINGULARITY_CACHEDIR'], ignore_errors=True)
        raise
    log_info('Built singularity sandbox for {} in {} seconds'.format(image, round(time.time() - start_time, 4)))

    shutil.rmtree(download_env['SINGULARITY_CACHEDIR'], ignore_errors=True)

    try:
        os.rename(temp_sandbox_dir, sandbox_dir)
    except OSError as e:
        if e.errno in (errno.EEXIST, errno.ENOTEMPTY):
            # somebody else won the race, which is fine -- drop our redundant copy
            assert os.path.exists(sandbox_dir)
            shutil.rmtree(temp_sandbox_dir, ignore_errors=True)
        else:
            raise

    return sandbox_dir


def prepare_image(image, work_dir=None, file_store=None):
    """ pull an image up front so that the fan-out does not have every job pulling at once """
    mode = binaries_mode()
    if mode == 'singularity':
        singularity_sandbox(image, work_dir=work_dir, file_store=file_store)
    elif mode == 'docker':
        log_info('Pulling docker image {}'.format(image))
        subprocess.check_call(['docker', 'pull', image])


def check_bind_mount(image, work_dir):
    """ make sure the container runtime can actually see the work directory.

    A snap-confined docker silently bind-mounts an *empty* directory when the host path is
    outside $HOME -- no error, no warning.  The first sign of trouble is a tool reporting that
    its inputs do not exist, several jobs into a fan-out, which is a miserable thing to debug.
    Probe for it once, up front, and say what to do about it. """
    mode = binaries_mode()
    if mode == 'local':
        return

    probe_dir = os.path.join(work_dir, 'bind-probe')
    os.makedirs(probe_dir, exist_ok=True)
    sentinel = 'toil-vcfcomp-bind-probe'
    with open(os.path.join(probe_dir, sentinel), 'w') as probe_file:
        probe_file.write('ok\n')

    # go through a shell so an image without its own `ls` is not mistaken for a mount failure,
    # and so a runtime that fails to start is reported as itself rather than excused
    try:
        listing = run_tool(['sh', '-c', 'ls ' + CONTAINER_MNT], image=image, work_dir=probe_dir,
                           check_output=True)
    except FileNotFoundError as e:
        log_info('Skipping the bind-mount check for {}: {}'.format(image, e))
        return
    except subprocess.CalledProcessError as e:
        raise RuntimeError(
            'Could not run a probe command in {} under {}: {}. The comparison jobs would all fail '
            'the same way.'.format(image, mode, e)) from e

    if sentinel not in listing:
        raise RuntimeError(
            'The {} runtime mounted {} at {} as an empty directory. Files written there are '
            'invisible inside the container, so every comparison would fail with "file does not '
            'exist". This is what a snap-confined docker does with any path outside $HOME. Rerun '
            'with --workDir pointing somewhere the container runtime can read (eg a directory '
            'under your home), or use a docker installed outside snap.'.format(
                mode, probe_dir, CONTAINER_MNT))


def check_host_tools(tools):
    """ fail early and with an actionable message rather than deep inside a fan-out """
    missing = [t for t in tools if not shutil.which(t)]
    if missing:
        raise RuntimeError('Required tool(s) not found on PATH: {}. toil-vcfcomp needs {} '
                           'available on every worker.'.format(', '.join(missing), ', '.join(tools)))
