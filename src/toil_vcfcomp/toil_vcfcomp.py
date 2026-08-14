#!/usr/bin/env python3

"""
Compare every sample of a multi-sample VCF against its own truth set, in parallel with Toil.

The query is typically a pangenome graph VCF with a few hundred samples in it; the truth is a
manifest (from hprc-truthset) naming one VCF and one confident-region BED per sample.  Samples
present in both are compared one at a time -- each comparison is an independent job -- and the
per-sample results are collected into a single report.

Several query VCFs can be given at once.  Each is a callset, benchmarked against the same truth
and reported side by side, so graphs can be compared directly rather than by diffing two runs.
Each sample's truth is fetched once and shared by every callset.  A callset is named after its
file, or explicitly with --queryVcf NAME=path, which is worth doing when the files are named after
their release directories and so collide.

Example:

  toil-vcfcomp ./jobstore \
      --queryVcf hprc-r2-mc-grch38.vcf.gz hprc-r2-mc-grch38.wave.vcf.gz \
      --truthIndex truth.index.csv \
      --ref GRCh38_no_alt.fa \
      --outDir ./report \
      --evalType aardvark \
      --batchSystem slurm --binariesMode singularity
"""

import os, sys
import subprocess
import timeit
import shutil
import shlex
import glob
import traceback
from urllib.parse import urlparse

from toil.common import Toil
from toil.job import Job
from toil.statsAndLogging import logger, set_logging_from_options
from toil.realtimeLogger import RealtimeLogger
from toil.lib.humanize import bytes2human, human2bytes

from toil_vcfcomp import vcfcomp
from toil_vcfcomp import report
from toil_vcfcomp.common import (makeURL, add_binaries_options, setup_binaries, prepare_image,
                                 check_host_tools, check_bind_mount, run_tool, binaries_mode, MNT)
from toil_vcfcomp.manifest import read_manifest, download_url

try:
    from toil_vcfcomp.version import toil_vcfcomp_commit
except ImportError:
    toil_vcfcomp_commit = 'unknown'

# the tools each write their headline numbers to a different set of files.  these are the ones
# worth keeping for every sample; --exportAll gets everything else
summary_globs = {'aardvark': ['*.tsv'],
                 'vcfeval': ['summary.txt', 'phasing.txt'],
                 'truvari': ['summary.json', 'summary.txt', 'params.json'],
                 'happy': ['*.summary.csv', '*.extended.csv']}

tool_dockers = {'aardvark': vcfcomp.default_aardvark_docker,
                'vcfeval': vcfcomp.default_vcfeval_docker,
                'truvari': vcfcomp.default_truvari_docker,
                'happy': vcfcomp.default_happy_docker}

# what each tool is called on the PATH, for --binariesMode local
local_binaries = {'aardvark': 'aardvark', 'vcfeval': 'rtg', 'truvari': 'truvari', 'happy': 'hap.py'}


def human2bytesN(s):
    if s is None:
        return None
    n = human2bytes(str(s))
    if n < 10000000:
        raise RuntimeError('Suspiciously small disk or memory specification: {}. '
                           'Did you forget a "Gi"?'.format(s))
    return n


def clamp_memory(memory_bytes):
    """ keep a computed memory request inside the --maxMemory / --defaultMemory the user gave """
    return max(min(int(os.environ['VCFCOMP_MAX_MEMORY']), int(memory_bytes)),
               int(os.environ['VCFCOMP_DEFAULT_MEMORY']))


def override_toil_options(options):
    """ toil defaults that make sense here, plus stashing the memory limits where clamp_memory()
    can reach them without carrying options around """
    if options.retryCount is None and options.batchSystem.lower() not in ['single_machine', 'singleMachine']:
        options.retryCount = 5

    if not options.realTimeLogging:
        logger.info('Enabling realtime logging in Toil')
        options.realTimeLogging = True

    max_mem = human2bytes(str(options.maxMemory)) if options.maxMemory else sys.maxsize
    if options.batchSystem.lower() in ['single_machine', 'singleMachine']:
        from toil import physicalMemory
        max_mem = min(max_mem, physicalMemory())
    os.environ['VCFCOMP_MAX_MEMORY'] = str(max_mem)
    os.environ['VCFCOMP_DEFAULT_MEMORY'] = str(
        human2bytes(str(options.defaultMemory)) if options.defaultMemory else 2 ** 31)


# --- jobs ---

def vcfcomp_workflow(job, options, manifest_rows, query_ids, ref_id, ref_fai_id):
    """ preflight, then split the query and build the SDF, then fan out over samples """
    # prepare_tools genuinely has to finish first: make_sdf runs rtg in a container, so it must
    # not start before the image is cached and before check_bind_mount has had its say -- otherwise
    # an unusable mount shows up as rtg's opaque "file does not exist" instead of our message
    prep_job = job.addChildJobFn(prepare_tools, options)
    root_job = Job()
    prep_job.addFollowOn(root_job)

    truth_samples = sorted(set(row['sample'] for row in manifest_rows))
    split_rvs = {}
    for callset, path, vcf_id, tbi_id in query_ids:
        split_job = root_job.addChildJobFn(split_query, options, callset, path, vcf_id, tbi_id,
                                           truth_samples,
                                           disk=split_disk(options, vcf_id, len(truth_samples)),
                                           memory=split_memory(options, len(truth_samples)))
        split_rvs[callset] = split_job.rv()

    # everything the reference needs -- decompressing, indexing, and the rtg SDF -- happens in one
    # job on a worker, so a remote or gzipped --ref never has to be prepared by hand first.  it
    # runs alongside the splits rather than after them.
    need_sdf = 'vcfeval' in options.evalType or 'happy' in options.evalType
    ref_job = root_job.addChildJobFn(prepare_reference, options, ref_id, ref_fai_id, need_sdf,
                                     disk=reference_disk(options, ref_id, need_sdf),
                                     memory=clamp_memory(8 * 2 ** 30))

    # .rv('fasta') rather than .rv(): by the time compare_all_samples runs these are real FileIDs,
    # so the per-comparison disk can be sized from the *uncompressed* reference
    fan_job = root_job.addFollowOnJobFn(compare_all_samples, options, manifest_rows, split_rvs,
                                        ref_job.rv('fasta'), ref_job.rv('fai'), ref_job.rv('sdf'))
    return fan_job.rv()


def split_memory(options, n_truth_samples):
    """ memory for one split job.

    bcftools +split scales with the number of samples it writes out, not with the input size.
    Measured on chr1 of a 232-sample pangenome VCF: splitting 231 samples peaked at 9.8 GB and took
    52 minutes, splitting one sample peaked at 0.35 GB and took 39 seconds.  So a run against a
    single-sample truth set (GIAB, say) needs a fraction of what a whole-cohort run does, and a
    flat default is wrong for one or the other.

    Scaling within a sample count is superlinear in the records handled and has not been measured
    at whole-genome scale, so this is a floor rather than a guarantee -- pair it with --doubleMem.
    """
    if options.splitMemory:
        return clamp_memory(options.splitMemory)
    return clamp_memory(max(4 * 2 ** 30, int(0.1 * 2 ** 30 * max(1, n_truth_samples))))


def split_disk(options, vcf_id, n_truth_samples):
    """ disk for one split job.

    Splitting a multi-sample VCF into single-sample ones is not size-preserving: every output
    repeats the site columns, so the total grows with the number of samples split out rather than
    staying near the input size.  Measured on a real 232-sample pangenome VCF the outputs came to
    11.4x the input for 231 samples -- about 0.05x the input per sample -- against the flat 4x
    this used to request, which was taken from a toy file and was wrong by 2.5x at real scale.
    """
    if options.splitDisk:
        return options.splitDisk
    per_sample = 0.05 * max(1, n_truth_samples)
    # 1x for the local copy of the input, the rest for the outputs
    return int(vcf_id.size * (1.2 + per_sample)) + 4 * 2 ** 30


def reference_disk(options, ref_id, need_sdf):
    """ disk for the reference job.  the imported size is all we know at graph-construction time,
    and for a gzipped fasta that is about a third of what gets unpacked """
    expanded = ref_id.size * 4 if options.ref.endswith('.gz') else ref_id.size
    total = ref_id.size + expanded * 1.2
    if need_sdf:
        # the SDF plus the tar of it
        total += expanded * 2.5
    return int(total) + 4 * 2 ** 30


def prepare_tools(job, options):
    """ preflight: fail now with a clear message rather than in the middle of a fan-out, and get
    the container images into the cache once instead of a few hundred times at once """
    check_host_tools(['bcftools', 'tabix', 'bgzip'])
    if options.regions:
        check_host_tools(['bedtools'])
    if any(row_url.startswith('s3://') for row_url in getattr(options, 'truth_urls', [])):
        # the s3:// form of a manifest is fetched with the AWS CLI, which is not a dependency of
        # this package and is routinely absent from compute nodes
        check_host_tools(['aws'])

    if binaries_mode() == 'local':
        # nothing containerises the comparators in local mode, so they have to be on the PATH
        check_host_tools(sorted(set(local_binaries[tool] for tool in options.evalType)))

    work_dir = job.fileStore.getLocalTempDir()
    for tool in options.evalType:
        prepare_image(options.docker[tool], work_dir=work_dir, file_store=job.fileStore)
        # probe every image, not just the first: they do not all have the same filesystem layout
        check_bind_mount(options.docker[tool], work_dir)


def prepare_reference(job, options, ref_id, ref_fai_id, need_sdf):
    """ get the reference into the shape every comparison needs, once, on a worker.

    --ref may be a URL and may be gzipped, so that a 3 GB fasta never has to be fetched, unpacked
    and indexed by hand before the run.  Decompressing, indexing and building the rtg SDF all
    happen here rather than in three jobs, because each step needs the output of the last on the
    same disk and there is nothing to gain from moving a multi-GB fasta between jobs to do it.
    """
    check_host_tools(['samtools'])
    work_dir = job.fileStore.getLocalTempDir()

    local = os.path.join(work_dir, os.path.basename(options.ref))
    job.fileStore.readGlobalFile(ref_id, local, mutable=True)

    if local.endswith('.gz'):
        RealtimeLogger.info('Decompressing reference {}'.format(os.path.basename(options.ref)))
        # bgzipped is fine for gzip -d too, and the comparison tools all need it plain
        subprocess.check_call(['gzip', '-fd', local])
        local = local[:-len('.gz')]
    ref_path = os.path.join(work_dir, options.refName)
    if os.path.abspath(local) != os.path.abspath(ref_path):
        os.rename(local, ref_path)

    if ref_fai_id is not None:
        job.fileStore.readGlobalFile(ref_fai_id, ref_path + '.fai')
    else:
        RealtimeLogger.info('Indexing reference {}'.format(options.refName))
        subprocess.check_call(['samtools', 'faidx', ref_path])

    check_reference_contigs(options, ref_path + '.fai')

    out = {'fasta': job.fileStore.writeGlobalFile(ref_path),
           'fai': job.fileStore.writeGlobalFile(ref_path + '.fai'),
           'sdf': None}

    if need_sdf:
        # built once here; every vcfeval/hap.py comparison would otherwise redo it
        RealtimeLogger.info('Building rtg SDF for {}'.format(options.refName))
        run_tool(['rtg', 'format', MNT + options.refName, '-o', MNT + 'ref.SDF'],
                 image=options.docker['vcfeval'], work_dir=work_dir, file_store=job.fileStore)
        sdf_tar = os.path.join(work_dir, 'ref.SDF.tar')
        subprocess.check_call(['tar', '-cf', sdf_tar, '-C', work_dir, 'ref.SDF'])
        out['sdf'] = job.fileStore.writeGlobalFile(sdf_tar)

    return out


def check_reference_contigs(options, fai_path):
    """ the reference-vs-query contig check, run where the .fai actually exists.

    bcftools norm -c w passes a reference-mismatched record through unchanged and exits 0, so the
    wrong --ref produces a completed run full of confident wrong numbers.  options.queryContigs is
    read from the query VCF headers on the leader, which works for a URL as well as a local file.
    """
    ref_contigs = {}
    with open(fai_path) as fai_file:
        for line in fai_file:
            fields = line.split('\t')
            if len(fields) >= 2:
                ref_contigs[fields[0]] = int(fields[1])

    for label, contigs in sorted(getattr(options, 'queryContigs', {}).items()):
        shared = set(contigs) & set(ref_contigs)
        if not shared:
            raise RuntimeError(
                '{} shares no contig names with --ref ({}). It has contigs like {}, the reference '
                'has {}. These are different assemblies.'.format(
                    label, options.refName, ', '.join(sorted(contigs)[:5]),
                    ', '.join(sorted(ref_contigs)[:5])))
        mismatched = sorted(c for c in shared if contigs[c] != ref_contigs[c])
        if mismatched:
            raise RuntimeError(
                '{} disagrees with --ref ({}) on the length of {} contig(s), eg {}. Comparing '
                'against the wrong reference produces plausible but wrong numbers, so this is '
                'fatal.'.format(label, options.refName, len(mismatched),
                                ', '.join('{} ({} vs {})'.format(c, contigs[c], ref_contigs[c])
                                          for c in mismatched[:3])))
        RealtimeLogger.info('{} agrees with the reference on {} shared contig(s)'.format(
            label, len(shared)))


def split_query(job, options, callset, query_vcf, query_vcf_id, query_tbi_id, truth_samples):
    """ one pass over the big query VCF, writing one single-sample VCF per shared sample.  doing
    this with a bcftools view per sample would mean reading the whole thing once per sample """
    work_dir = job.fileStore.getLocalTempDir()
    query_path = os.path.join(work_dir, os.path.basename(query_vcf))
    job.fileStore.readGlobalFile(query_vcf_id, query_path)
    if query_tbi_id:
        job.fileStore.readGlobalFile(query_tbi_id, query_path + '.tbi')

    query_samples = subprocess.check_output(['bcftools', 'query', '-l', query_path],
                                            encoding='utf-8').split()
    shared = sorted(set(query_samples) & set(truth_samples))
    if options.samples:
        shared = sorted(set(shared) & set(options.samples))

    query_only = sorted(set(query_samples) - set(truth_samples))
    truth_only = sorted(set(truth_samples) - set(query_samples))
    RealtimeLogger.info('{}: {} samples shared between the query VCF ({}) and the truth index ({})'.format(
        callset, len(shared), len(query_samples), len(truth_samples)))
    # a sample missing from one side is nearly always a name mismatch rather than a real absence,
    # so name them instead of silently comparing fewer samples than the user expected
    if query_only:
        RealtimeLogger.warning('{}: {} sample(s) only in the query VCF, not compared: {}'.format(
            callset, len(query_only), ', '.join(query_only[:10]) + (' ...' if len(query_only) > 10 else '')))
    if truth_only:
        RealtimeLogger.warning('{}: {} sample(s) only in the truth index, not compared: {}'.format(
            callset, len(truth_only), ', '.join(truth_only[:10]) + (' ...' if len(truth_only) > 10 else '')))
    if not shared:
        raise RuntimeError('No samples in common between {} and {}. The query VCF has samples like '
                           '{}, and the truth index has samples like {}.'.format(
                               query_vcf, options.truthIndex,
                               ', '.join(sorted(query_samples)[:5]), ', '.join(truth_samples[:5])))

    samples_path = os.path.join(work_dir, 'shared-samples.txt')
    with open(samples_path, 'w') as samples_file:
        for sample in shared:
            samples_file.write(sample + '\n')

    split_dir = os.path.join(work_dir, 'split')
    os.makedirs(split_dir)
    split_cmd = ['bcftools', '+split', query_path, '-S', samples_path, '-Oz', '-o', split_dir,
                 '-i', options.splitFilter]
    if options.splitKeepTags:
        split_cmd += ['-k', options.splitKeepTags]
    RealtimeLogger.info('Splitting {} into {} single-sample VCFs'.format(query_vcf, len(shared)))
    subprocess.check_call(split_cmd)

    sample_ids = {}
    for sample in shared:
        sample_vcf = os.path.join(split_dir, sample + '.vcf.gz')
        if not os.path.isfile(sample_vcf):
            RealtimeLogger.warning('bcftools +split produced no output for {} -- it has no non-ref '
                                   'genotypes in {}, so it is skipped'.format(sample, callset))
            continue
        sample_ids[sample] = job.fileStore.writeGlobalFile(sample_vcf)
    return sample_ids


def compare_all_samples(job, options, manifest_rows, split_by_callset, ref_id, ref_fai_id, sdf_id):
    """ fan out over callset x sample x tool.

    Truth is fetched once per sample and reused by every callset: the callsets are being compared
    against the same truth, so downloading it once per callset would multiply the transfer for
    nothing.
    """
    # union of the samples each callset actually has, so callsets need not agree on their sample sets
    all_samples = set(s for ids in split_by_callset.values() for s in ids)
    # one truth row per (sample, truthset): a sample can be benchmarked against several truth sets
    truth_rows = [row for row in manifest_rows if row['sample'] in all_samples]
    truthsets = sorted(set(row['truthset'] for row in truth_rows))
    RealtimeLogger.info('Comparing {} callset(s) against {} truthset(s) ({}) over {} sample(s) '
                        'with {}'.format(len(split_by_callset), len(truthsets), ', '.join(truthsets),
                                         len(all_samples), ', '.join(options.evalType)))
    all_samples = sorted(all_samples)
    for callset in sorted(split_by_callset):
        missing = sorted(set(all_samples) - set(split_by_callset[callset]))
        if missing:
            RealtimeLogger.warning('callset {} is missing {} of the {} samples, so its rows will '
                                   'cover fewer samples than the others: {}'.format(
                                       callset, len(missing), len(all_samples),
                                       ', '.join(missing[:10]) + (' ...' if len(missing) > 10 else '')))

    results = []
    # one fetch per (sample, truthset) -- shared by every callset, since they are all scored
    # against the same truth
    for row in sorted(truth_rows, key=lambda r: (r['sample'], r['truthset'])):
        sample, truthset = row['sample'], row['truthset']
        fetch_job = job.addChildJobFn(fetch_truth, options, row,
                                      disk=4 * 2 ** 30, memory=2 ** 30, cores=1)
        for callset in sorted(split_by_callset):
            if sample not in split_by_callset[callset]:
                continue
            for tool in options.evalType:
                compare_job = fetch_job.addFollowOnJobFn(compare_sample, options, callset, truthset,
                                                         sample, tool,
                                                         split_by_callset[callset][sample], fetch_job.rv(),
                                                         ref_id, ref_fai_id, sdf_id,
                                                         cores=options.threads,
                                                         memory=clamp_memory(options.compareMemory),
                                                         disk=int(ref_id.size * 5) + options.compareDisk)
                results.append(compare_job.rv())

    report_job = job.addFollowOnJobFn(make_report, options, results)
    return report_job.rv()


def fetch_truth(job, options, row):
    """ get one sample's truth VCF, index and BED into the job store.  done in its own small job
    so a few hundred downloads run in parallel instead of serially through the leader """
    work_dir = job.fileStore.getLocalTempDir()

    def stage(source):
        local = os.path.join(work_dir, os.path.basename(source))
        if os.path.isfile(source):
            return source
        RealtimeLogger.info('Fetching {}'.format(source))
        download_url(source, local)
        return local

    if not row['vcf']:
        raise RuntimeError('truth index row for {} / {} has no VCF'.format(
            row['sample'], row.get('truthset')))

    ids = {}
    vcf_local = stage(row['vcf'])
    ids['vcf'] = job.fileStore.writeGlobalFile(vcf_local)

    # not every published benchmark ships a .tbi, and the workflow needs one to read contigs and
    # infer sex, so build it rather than failing
    try:
        tbi_local = stage(row['vcf'] + '.tbi')
    except Exception as e:
        RealtimeLogger.info('No .tbi for {} ({}), indexing it here instead'.format(row['vcf'], e))
        if vcf_local == row['vcf']:
            # the source is a local file we must not write next to; index a copy
            vcf_local = os.path.join(work_dir, os.path.basename(row['vcf']))
            shutil.copyfile(row['vcf'], vcf_local)
        subprocess.check_call(['tabix', '-fp', 'vcf', vcf_local])
        tbi_local = vcf_local + '.tbi'
    ids['vcf.tbi'] = job.fileStore.writeGlobalFile(tbi_local)

    if row['bed']:
        bed_local = stage(row['bed'])
        if bed_local.endswith('.gz'):
            # some GIAB benchmarks ship the confident regions gzipped, and every comparison tool
            # wants a plain BED
            if bed_local == row['bed']:
                copied = os.path.join(work_dir, os.path.basename(row['bed']))
                shutil.copyfile(row['bed'], copied)
                bed_local = copied
            subprocess.check_call(['gzip', '-fd', bed_local])
            bed_local = bed_local[:-len('.gz')]
        ids['bed'] = job.fileStore.writeGlobalFile(bed_local)

    return ids


def infer_haploid_sex(vcf_path):
    """ chrY records are present if and only if the sample is male, so the truth VCF tells us
    whether chrX should be treated as haploid without needing any sample metadata """
    contigs = subprocess.check_output(['bcftools', 'index', '-s', vcf_path], encoding='utf-8')
    with_records = set(line.split('\t')[0].lower() for line in contigs.strip().split('\n') if line.strip())
    return 'chry' in with_records


def compare_sample(job, options, callset, truthset, sample, tool, query_vcf_id, truth_ids,
                   ref_id, ref_fai_id, sdf_id):
    """ the actual comparison for one sample with one tool.

    Never propagates a failure.  make_report is a follow-on of the whole fan-out, so a single job
    that exhausts its retries would take the entire report down with it -- after every other
    sample had already been paid for.  Failures are recorded and returned instead, and main()
    decides at the end whether there were too many.
    """
    try:
        return run_comparison(job, options, callset, truthset, sample, tool, query_vcf_id,
                              truth_ids, ref_id, ref_fai_id, sdf_id)
    except Exception as e:
        RealtimeLogger.error('{}/{}/{}/{} FAILED: {}\n{}'.format(callset, truthset, sample, tool, e,
                                                                 traceback.format_exc()))
        return {'rows': [], 'chrom_rows': [], 'callset': callset, 'truthset': truthset,
                'sample': sample, 'tool': tool, 'error': str(e)}


def run_comparison(job, options, callset, truthset, sample, tool, query_vcf_id, truth_ids,
                   ref_id, ref_fai_id, sdf_id):
    work_dir = job.fileStore.getLocalTempDir()
    out_dir = os.path.join(work_dir, tool)
    os.makedirs(out_dir)

    ref_path = os.path.join(work_dir, options.refName)
    job.fileStore.readGlobalFile(ref_id, ref_path)
    job.fileStore.readGlobalFile(ref_fai_id, ref_path + '.fai')

    truth_vcf = os.path.join(work_dir, '{}.truth.vcf.gz'.format(sample))
    job.fileStore.readGlobalFile(truth_ids['vcf'], truth_vcf)
    job.fileStore.readGlobalFile(truth_ids['vcf.tbi'], truth_vcf + '.tbi')

    query_vcf = os.path.join(work_dir, '{}.query.vcf.gz'.format(sample))
    job.fileStore.readGlobalFile(query_vcf_id, query_vcf)
    subprocess.check_call(['tabix', '-fp', 'vcf', query_vcf])

    regions_bed = resolve_regions(job, options, sample, truth_ids, work_dir)

    male = infer_haploid_sex(truth_vcf)
    haploid_x = male
    remove_y = not male
    RealtimeLogger.info('{}: chrY {} in the truth set, so running with haploid_x={} remove_y={}'.format(
        sample, 'present' if male else 'absent', haploid_x, remove_y))

    min_length, max_length, bi_allelic = preprocess_params(options, tool)

    # unlike the GIAB truth sets vcfcomp was written against, dipcall output is not preconditioned:
    # it is unnormalized, and hemizygous chrX/chrY come through as .|1 / 1|. .  so both sides get
    # the same treatment here
    par_regions = vcfcomp.parse_par_bed(options.parBed) if options.parBed else None

    truth_prepped = os.path.join(out_dir, '{}.truth.prep.vcf.gz'.format(sample))
    vcfcomp.vcf_preprocess(truth_vcf, truth_prepped, ref_path, remove_y, haploid_x, sample,
                           min_length, max_length, bi_allelic, True,
                           filter_expr=options.truthFilter, par_regions=par_regions)

    query_prepped = os.path.join(out_dir, '{}.query.prep.vcf.gz'.format(sample))
    vcfcomp.vcf_preprocess(query_vcf, query_prepped, ref_path, remove_y, haploid_x, sample,
                           min_length, max_length, bi_allelic, True,
                           filter_expr=options.queryFilter, par_regions=par_regions)

    sdf_path = untar_sdf(job, sdf_id, work_dir, out_dir) if sdf_id else None

    if tool == 'aardvark':
        vcfcomp.aardvark(truth_prepped, query_prepped, ref_path, regions_bed, sample, out_dir,
                         threads=int(job.cores), options=options.toolOptions[tool],
                         docker_image=options.docker[tool])
        chrom_table = vcfcomp.aardvark_chromsplit(out_dir)
    elif tool == 'vcfeval':
        vcfcomp.vcfeval(truth_prepped, query_prepped, ref_path, regions_bed, sample, out_dir,
                        threads=int(job.cores), options=options.toolOptions[tool],
                        docker_image=options.docker[tool], sdf=sdf_path)
        chrom_table = vcfcomp.vcfeval_chromsplit(out_dir)
    elif tool == 'truvari':
        vcfcomp.truvari(truth_prepped, query_prepped, ref_path, regions_bed, sample, out_dir,
                        options=options.toolOptions[tool], docker_image=options.docker[tool])
        chrom_table = vcfcomp.truvari_chromsplit(out_dir)
    elif tool == 'happy':
        vcfcomp.happy(truth_prepped, query_prepped, ref_path, regions_bed, sample, out_dir,
                      threads=int(job.cores), options=set_happy_gender(options.toolOptions[tool], male),
                      docker_image=options.docker[tool], sdf=sdf_path)
        chrom_table = vcfcomp.happy_chromsplit(os.path.join(out_dir, 'happy-output.vcf.gz'))
    else:
        raise RuntimeError('unknown --evalType: {}'.format(tool))

    rows = report.parse_summary(tool, out_dir, callset, truthset, sample)
    if not rows:
        RealtimeLogger.warning('{}/{}/{}/{}: no summary file found to parse -- the comparison ran '
                               'but produced nothing to report'.format(callset, truthset, sample, tool))
    chrom_rows = report.chromsplit_to_long(chrom_table, callset, truthset, sample, tool)

    export_sample_outputs(job, options, callset, truthset, sample, tool, out_dir)

    return {'rows': rows, 'chrom_rows': chrom_rows}


def untar_sdf(job, sdf_id, work_dir, out_dir):
    """ unpack the shared SDF straight into out_dir, so the tool wrappers find it already staged
    and do not spend another few GB copying it """
    sdf_tar = os.path.join(work_dir, 'ref.SDF.tar')
    job.fileStore.readGlobalFile(sdf_id, sdf_tar)
    subprocess.check_call(['tar', '-xf', sdf_tar, '-C', out_dir])
    os.remove(sdf_tar)
    return os.path.join(out_dir, 'ref.SDF')


def set_happy_gender(happy_options, male):
    """ hap.py's --gender changes how it treats the sex chromosomes.  a fixed value is fine when
    benchmarking one sample, but across a cohort it has to follow each sample """
    gender = 'male' if male else 'female'
    fields = [f for f in happy_options.split() if not f.startswith('--gender')]
    return ' '.join(fields + ['--gender={}'.format(gender)])


def preprocess_params(options, tool):
    """ mirror the per-tool preprocessing that vcfcomp's own subcommands use """
    if tool == 'happy':
        return 0, options.maxLength, True
    if tool == 'truvari':
        return options.minLength, options.svMaxLength, True
    if tool == 'vcfeval':
        return 0, options.maxLength, False
    return options.minLength, options.svMaxLength, False


def resolve_regions(job, options, sample, truth_ids, work_dir):
    """ the sample's own confident regions, optionally cut down by a global BED """
    bed_path = None
    if 'bed' in truth_ids:
        bed_path = os.path.join(work_dir, '{}.truth.bed'.format(sample))
        job.fileStore.readGlobalFile(truth_ids['bed'], bed_path)

    global_bed = None
    if options.regions:
        global_bed = os.path.join(work_dir, 'regions.bed')
        job.fileStore.readGlobalFile(options.regions_id, global_bed)

    if bed_path and global_bed:
        intersected = os.path.join(work_dir, '{}.regions.bed'.format(sample))
        # bedtools intersect emits one interval per overlapping (a,b) pair, so a --regions BED that
        # overlaps itself -- which repeat and segdup annotations routinely do -- yields duplicated
        # evaluation intervals.  merge, or the doubled spans get double-counted in the denominators
        with open(intersected, 'w') as bed_file:
            vcfcomp.run_pipe('bedtools intersect -a {} -b {} | sort -k1,1 -k2,2n | bedtools merge'.format(
                shlex.quote(bed_path), shlex.quote(global_bed)), stdout=bed_file)
        regions = intersected
    elif bed_path:
        regions = bed_path
    elif global_bed:
        regions = global_bed
    else:
        raise RuntimeError('sample {} has no confident-region BED in the truth index, and no '
                           '--regions was given. Every comparison tool needs one.'.format(sample))

    if os.path.getsize(regions) == 0:
        raise RuntimeError('sample {} has an empty evaluation region set. If --regions was given, '
                           'it does not overlap that sample\'s confident regions at all.'.format(sample))
    return regions


def export_sample_outputs(job, options, callset, truthset, sample, tool, out_dir):
    """ export the per-sample summaries (and, with --exportAll, the whole output directory) """
    if not options.outDir:
        return
    dest_dir = os.path.join(options.outDir, callset, truthset, sample, tool)
    if not options.outDir.startswith('s3://') and not os.path.isdir(dest_dir):
        os.makedirs(dest_dir, exist_ok=True)

    for pattern in summary_globs.get(tool, []):
        for path in sorted(glob.glob(os.path.join(out_dir, pattern))):
            file_id = job.fileStore.writeGlobalFile(path)
            job.fileStore.exportFile(file_id, makeURL(os.path.join(dest_dir, os.path.basename(path))))
            # the job store copy exists only to hand the file to exportFile
            job.fileStore.deleteGlobalFile(file_id)

    if options.exportAll:
        # the reference is hard-linked into out_dir so the container can see it, which costs no
        # disk but which tar reads in full regardless of link count -- on a whole human genome
        # that is ~94% of the archive, repeated for every sample.  the SDF is the same story.
        sdf_dir = os.path.join(out_dir, 'ref.SDF')
        if os.path.isdir(sdf_dir):
            shutil.rmtree(sdf_dir)
        ref_name = options.refName
        tar_path = os.path.join(os.path.dirname(out_dir),
                                '{}.{}.{}.{}.tar'.format(callset, truthset, sample, tool))
        # no -z: what is left is bgzf-compressed already, so gzip would burn CPU for nothing
        subprocess.check_call(['tar', '-cf', tar_path,
                               '--exclude', ref_name, '--exclude', ref_name + '.fai',
                               '-C', os.path.dirname(out_dir), os.path.basename(out_dir)])
        file_id = job.fileStore.writeGlobalFile(tar_path)
        job.fileStore.exportFile(file_id, makeURL(os.path.join(dest_dir, os.path.basename(tar_path))))
        job.fileStore.deleteGlobalFile(file_id)


def make_report(job, options, results):
    """ collect every sample's rows into the three report tables """
    work_dir = job.fileStore.getLocalTempDir()

    rows = []
    chrom_rows = []
    failures = []
    for result in results:
        if not result:
            continue
        if result.get('error'):
            failures.append('{}/{}/{}/{}: {}'.format(result.get('callset'), result.get('truthset'),
                                                     result.get('sample'), result.get('tool'),
                                                     result['error']))
            continue
        rows += result['rows']
        chrom_rows += result['chrom_rows']

    if failures:
        RealtimeLogger.warning('{} comparison(s) failed and are absent from the report:\n  {}'.format(
            len(failures), '\n  '.join(failures)))

    if not rows:
        raise RuntimeError('No results to report: every comparison failed or produced an '
                           'unparseable summary. Check the job logs for the individual samples.')

    paths = report.write_report(rows, chrom_rows, work_dir,
                                callset_order=getattr(options, 'callsets', None))
    RealtimeLogger.info('Report covers {} rows from {} callset/truthset/sample/tool comparisons'.format(
        len(rows), len(set((r['callset'], r['truthset'], r['sample'], r['tool']) for r in rows))))
    file_ids = dict((name, job.fileStore.writeGlobalFile(path)) for name, path in paths.items())
    file_ids['_failures'] = failures
    return file_ids


# --- main ---

def vcf_header_contigs(vcf_path):
    """ {contig: length} from a VCF header.  works on a URL as well as a local file """
    header = subprocess.check_output(['bcftools', 'view', '-h', vcf_path], encoding='utf-8')
    contigs = {}
    for line in header.split('\n'):
        if not line.startswith('##contig='):
            continue
        fields = dict(f.split('=', 1) for f in line[len('##contig=<'):].rstrip('>').split(',') if '=' in f)
        if 'ID' in fields and 'length' in fields:
            contigs[fields['ID']] = int(fields['length'])
    return contigs


# the same assembly under the names UCSC and the benchmark producers each use.  only used to
# decide whether to warn; the contig-length check runs regardless and is what actually catches a
# wrong reference
assembly_aliases = [set(['hs1', 'chm13v2', 'chm13v2.0']),
                    set(['hg38', 'grch38', 'grch38_no_alt'])]


def same_assembly(manifest_ref, ref_name):
    a = manifest_ref.strip().lower()
    b = os.path.splitext(ref_name)[0].strip().lower()
    if a in b or b in a:
        return True
    return any(a in group and b in group for group in assembly_aliases)


def collect_query_contigs(options, manifest_rows):
    """ read the contig lists the reference will be checked against.

    Only VCF *headers* are read, which is a few hundred KB even for a remote multi-GB VCF, so this
    stays cheap when --queryVcf is a URL.  The comparison itself happens in prepare_reference,
    where the .fai exists.
    """
    manifest_refs = set(row['ref'] for row in manifest_rows if row.get('ref'))
    for manifest_ref in sorted(manifest_refs):
        if manifest_ref and not same_assembly(manifest_ref, options.refName):
            logger.warning('The truth index says its calls are against "{}", but --ref is "{}". '
                           'Check that these are the same assembly.'.format(manifest_ref, options.refName))

    to_check = [('query VCF {}'.format(os.path.basename(p)), p) for p in options.queryVcf]
    # one truth VCF is enough to catch a wrong reference, and truth VCFs are often remote
    if manifest_rows:
        row = manifest_rows[0]
        to_check.append(('truth VCF {} / {}'.format(row['sample'], row['truthset']), row['vcf']))

    contigs = {}
    for label, vcf_path in to_check:
        try:
            contigs[label] = vcf_header_contigs(vcf_path)
        except Exception as e:
            logger.warning('Could not read the header of {} for the contig check: {}'.format(label, e))
    return contigs


def split_query_spec(spec):
    """ 'NAME=path' -> ('NAME', 'path'), anything else -> (None, spec).

    Only treated as a label when what precedes the first '=' looks like a name rather than part
    of a path, so a filename or URL containing '=' is left alone.
    """
    if '=' in spec:
        head, tail = spec.split('=', 1)
        if head and tail and '/' not in head and ':' not in head:
            return head, tail
    return None, spec


def callset_labels(specs):
    """ name each query VCF for the report.

    An explicit NAME=path wins.  Otherwise the basename without its extensions is the obvious
    label, but several callsets can share one (release directories tend to name the file after the
    directory), so a clash falls back to including the parent directory -- which is unambiguous but
    long, and a sign that NAME=path would read better.
    """
    def strip(path):
        name = os.path.basename(path)
        for suffix in ['.vcf.gz', '.vcf.bgz', '.bcf', '.vcf']:
            if name.endswith(suffix):
                return name[:-len(suffix)]
        return name

    given = [split_query_spec(s)[0] for s in specs]
    paths = [split_query_spec(s)[1] for s in specs]

    auto = [strip(p) for p in paths]
    out = []
    for explicit, label, path in zip(given, auto, paths):
        if explicit:
            out.append(explicit)
            continue
        if auto.count(label) > 1:
            parent = os.path.basename(os.path.dirname(os.path.abspath(path)))
            label = '{}.{}'.format(parent, label) if parent else label
        out.append(label)

    duplicates = sorted(set(l for l in out if out.count(l) > 1))
    if duplicates:
        raise SystemExit('Could not derive distinct callset names: {} appears more than once. '
                         'Name them explicitly with --queryVcf NAME=path.'.format(
                             ', '.join(duplicates)))
    return out


def callset_paths(specs):
    """ the paths, with any NAME= prefix removed """
    return [split_query_spec(s)[1] for s in specs]


def parse_sample_list(value):
    if not value:
        return None
    if os.path.isfile(value):
        with open(value) as sample_file:
            return sorted(set(line.strip() for line in sample_file if line.strip()))
    return sorted(set(s.strip() for s in value.split(',') if s.strip()))


def main():
    parser = Job.Runner.getDefaultArgumentParser()
    # Toil gives --retryCount an argparse default of 1, so the "only if the user did not say"
    # guard in override_toil_options could never fire and cluster runs got 1 retry, not 5
    parser.set_defaults(retryCount=None)

    parser.add_argument('--queryVcf', required=True, nargs='+',
                        help='multi-sample VCF(s) to evaluate (eg from vg deconstruct). give several '
                        'to benchmark them side by side against the same truth: each becomes a '
                        'callset in the report, named after its file')
    parser.add_argument('--truthIndex', required=True,
                        help='truth manifest CSV from hprc-truthset (sample,ref,caller,vcf,bed)')
    parser.add_argument('--ref', required=True,
                        help='reference fasta the VCFs were called against. may be a local path or a '
                        'URL, and may be gzipped: it is fetched, unpacked and indexed on a worker, '
                        'so nothing has to be prepared by hand first')
    parser.add_argument('--outDir', required=True,
                        help='output directory (or s3:// prefix) for the report and per-sample summaries')
    parser.add_argument('--evalType', nargs='+', default=['aardvark'],
                        choices=['aardvark', 'vcfeval', 'truvari', 'happy'],
                        help='comparison tool(s) to run, multiple allowed [aardvark]')
    parser.add_argument('--samples',
                        help='restrict to these samples (comma-separated, or a file with one per line)')
    parser.add_argument('--regions',
                        help='global BED to intersect with each sample confident-region BED')
    parser.add_argument('--threads', type=int, default=8,
                        help='threads (and cores requested) per comparison [8]')
    parser.add_argument('--minLength', type=int, default=vcfcomp.default_aardvark_min_length,
                        help='ignore sites whose alleles are all shorter than this [{}]'.format(
                            vcfcomp.default_aardvark_min_length))
    parser.add_argument('--maxLength', type=int, default=vcfcomp.default_happy_max_length,
                        help='max allele length for the small-variant tools (happy, vcfeval) [{}]'.format(
                            vcfcomp.default_happy_max_length))
    parser.add_argument('--svMaxLength', type=int, default=vcfcomp.default_aardvark_max_length,
                        help='max allele length for the SV-aware tools (aardvark, truvari) [{}]'.format(
                            vcfcomp.default_aardvark_max_length))
    parser.add_argument('--truthFilter', default=vcfcomp.default_filter_expr,
                        help='bcftools expression for the truth records to keep. the default drops '
                        'dipcall\'s HET/GAP flagged calls while keeping the unflagged ones, which are '
                        'FILTER="." rather than PASS ["{}"]'.format(vcfcomp.default_filter_expr))
    parser.add_argument('--queryFilter', default=None,
                        help='bcftools expression for the query records to keep [keep everything]')
    parser.add_argument('--splitFilter', default='GT~"[1-9]"',
                        help='bcftools expression for which records to keep when splitting the query '
                        'into single-sample VCFs. the default keeps any genotype carrying a non-ref '
                        'allele. do NOT use GT="alt" here: bcftools classifies a half-called genotype '
                        'like .|1 as "mis" rather than "alt", so it silently drops every hemizygous '
                        'call -- all of chrY and non-PAR chrX in males, and in a graph VCF around 40%% '
                        'of the calls overall [GT~"[1-9]"]')
    parser.add_argument('--splitKeepTags', default='FMT/GT',
                        help='tags bcftools +split keeps when splitting the query. the default drops '
                        'INFO, which in a deconstructed graph VCF is mostly the AT tag and often over '
                        'half the file, and which none of the comparison tools read [FMT/GT]')
    parser.add_argument('--parBed',
                        help='BED of pseudoautosomal regions, where a male is diploid and chrX must '
                        'not be flattened. defaults to the GRCh38 PARs; pass CHM13 ones here')
    parser.add_argument('--maxFailedSamples', type=int, default=0,
                        help='how many per-sample comparisons may fail before the run is called a '
                        'failure. the report is written either way, covering whatever succeeded [0]')
    parser.add_argument('--exportAll', action='store_true',
                        help='export each comparison\'s full output directory as a tarball, not just '
                        'its summary files')
    parser.add_argument('--aardvarkOptions', default=vcfcomp.default_aardvark_options,
                        help='options for aardvark (surround in quotes)')
    parser.add_argument('--vcfevalOptions', default=vcfcomp.default_vcfeval_options,
                        help='options for rtg vcfeval (surround in quotes)')
    parser.add_argument('--truvariOptions', default=vcfcomp.default_truvari_options,
                        help='options for truvari (surround in quotes)')
    parser.add_argument('--happyOptions', default=vcfcomp.default_happy_options,
                        help='options for hap.py (surround in quotes)')
    parser.add_argument('--aardvarkDocker', default=vcfcomp.default_aardvark_docker)
    parser.add_argument('--vcfevalDocker', default=vcfcomp.default_vcfeval_docker)
    parser.add_argument('--truvariDocker', default=vcfcomp.default_truvari_docker)
    parser.add_argument('--happyDocker', default=vcfcomp.default_happy_docker)
    parser.add_argument('--compareMemory', type=human2bytesN, default='32Gi',
                        help='memory for each comparison job [32Gi]')
    parser.add_argument('--compareDisk', type=human2bytesN, default='20Gi',
                        help='disk for each comparison job, on top of 5x the reference size [20Gi]')
    parser.add_argument('--splitMemory', type=human2bytesN, default=None,
                        help='memory for each query-split job. the default scales with the number of '
                        'truth samples (0.1Gi each, floor 4Gi), because bcftools +split scales with '
                        'the samples it writes out: on chr1 of a 232-sample pangenome VCF, 231 '
                        'samples peaked at 9.8Gi and one sample at 0.35Gi. pair with --doubleMem')
    parser.add_argument('--splitDisk', type=human2bytesN, default=None,
                        help='disk for each query-split job. the default scales with the number of '
                        'truth samples, because splitting one multi-sample VCF into N single-sample '
                        'ones repeats the site columns N times -- measured at 11x the input for 231 '
                        'samples, not a fraction of it')
    add_binaries_options(parser)

    options = parser.parse_args()

    setup_binaries(options)
    set_logging_from_options(options)

    # the comparison tools need an uncompressed, indexed fasta, but the user should not have to
    # produce one: a gzipped and/or remote --ref is unpacked and indexed on a worker
    options.refName = os.path.basename(urlparse(options.ref).path or options.ref)
    if options.refName.endswith('.gz'):
        options.refName = options.refName[:-len('.gz')]
    if not options.refName:
        raise SystemExit('could not work out a filename from --ref {}'.format(options.ref))

    if not options.outDir.startswith('s3://'):
        if not os.path.isdir(options.outDir):
            os.makedirs(options.outDir)
        options.outDir = os.path.abspath(options.outDir)

    if 'happy' in options.evalType and '--pass-only' in options.happyOptions:
        # dipcall leaves its good calls at FILTER="." rather than PASS, so --pass-only throws
        # essentially the whole truth set away and reports a recall near zero
        logger.warning('--happyOptions contains --pass-only. If the truth set is dipcall output, that '
                       'discards every record whose FILTER is "." -- which is most of them. Consider '
                       'dropping it and relying on --truthFilter instead.')

    # --queryVcf entries may be NAME=path; resolve both halves once, here
    options.callsets = callset_labels(options.queryVcf)
    options.queryVcf = callset_paths(options.queryVcf)

    options.samples = parse_sample_list(options.samples)
    options.docker = {'aardvark': options.aardvarkDocker, 'vcfeval': options.vcfevalDocker,
                      'truvari': options.truvariDocker, 'happy': options.happyDocker}
    options.toolOptions = {'aardvark': options.aardvarkOptions, 'vcfeval': options.vcfevalOptions,
                           'truvari': options.truvariOptions, 'happy': options.happyOptions}

    override_toil_options(options)

    # a request bigger than the batch system's ceiling is only rejected when the job is issued,
    # which for the comparison jobs is after the query VCF has been split and every truth file
    # downloaded.  catch it in the first second instead
    max_mem = int(os.environ['VCFCOMP_MAX_MEMORY'])
    for flag, value in [('--compareMemory', options.compareMemory), ('--splitMemory', options.splitMemory)]:
        # splitMemory defaults to None, meaning "derive it from the truth sample count"
        if value is not None and value > max_mem:
            raise SystemExit('{} is {} but the largest job this batch system will accept is {}. '
                             'Lower {} or raise --maxMemory.'.format(
                                 flag, bytes2human(value), bytes2human(max_mem), flag))

    logger.info('toil-vcfcomp Command: {}'.format(' '.join(sys.argv)))
    logger.info('toil-vcfcomp Commit: {}'.format(toil_vcfcomp_commit))
    start_time = timeit.default_timer()

    with Toil(options) as toil:
        if options.restart:
            report_id_dict = toil.restart()
        else:
            manifest_rows = read_manifest(options.truthIndex)
            logger.info('Read {} truth samples from {}'.format(len(manifest_rows), options.truthIndex))
            options.truth_urls = [row['vcf'] for row in manifest_rows]
            # headers only, so this is cheap even when everything is a URL.  the reference is
            # checked against these in prepare_reference, once the .fai exists
            options.queryContigs = collect_query_contigs(options, manifest_rows)

            query_ids = []
            for callset, query_vcf in zip(options.callsets, options.queryVcf):
                query_vcf_id = toil.importFile(makeURL(query_vcf))
                query_tbi_id = None
                if os.path.isfile(query_vcf + '.tbi'):
                    query_tbi_id = toil.importFile(makeURL(query_vcf + '.tbi'))
                query_ids.append((callset, query_vcf, query_vcf_id, query_tbi_id))
                logger.info('Callset "{}" is {} ({})'.format(
                    callset, os.path.basename(query_vcf), bytes2human(query_vcf_id.size)))
            ref_id = toil.importFile(makeURL(options.ref))
            ref_fai_id = None
            # only reuse an existing .fai when the reference is a plain local fasta; anything else
            # gets indexed on a worker, which is the whole point of allowing a URL here
            if not options.ref.endswith('.gz') and os.path.isfile(options.ref + '.fai'):
                ref_fai_id = toil.importFile(makeURL(options.ref + '.fai'))
            options.regions_id = None
            if options.regions:
                options.regions_id = toil.importFile(makeURL(options.regions))

            logger.info('Reference {} is {}{}'.format(
                options.refName, bytes2human(ref_id.size),
                ' compressed, and will be unpacked and indexed on a worker' if options.ref.endswith('.gz')
                else '' if ref_fai_id else ', and will be indexed on a worker'))

            report_id_dict = toil.start(Job.wrapJobFn(vcfcomp_workflow, options, manifest_rows,
                                                      query_ids, ref_id, ref_fai_id))

        failures = report_id_dict.pop('_failures', [])
        for name in sorted(report_id_dict):
            suffix = {'report': '.tsv', 'summary': '.summary.tsv', 'chrom': '.chrom.tsv',
                      'text': '.txt', 'compare': '.compare.tsv'}[name]
            out_path = makeURL(os.path.join(options.outDir, 'report' + suffix))
            logger.info('Exporting {}'.format(out_path))
            toil.exportFile(report_id_dict[name], out_path)

    end_time = timeit.default_timer()
    run_time = end_time - start_time
    logger.info('toil-vcfcomp has finished after {} seconds'.format(run_time))

    # decided only after the tables have been written and exported, so a handful of bad samples
    # never costs the report for all the good ones
    if failures:
        logger.warning('{} comparison(s) failed:\n  {}'.format(len(failures), '\n  '.join(failures)))
        if len(failures) > options.maxFailedSamples:
            raise SystemExit(
                '{} comparison(s) failed, more than --maxFailedSamples {}. The report covers the '
                'samples that succeeded and has been written to {}.'.format(
                    len(failures), options.maxFailedSamples, options.outDir))


if __name__ == '__main__':
    main()
