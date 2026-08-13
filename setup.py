from setuptools import setup, find_packages
import os
import subprocess

# stamp the git commit into the package so every run can log exactly what it was.  written at
# build time and gitignored, the way cactus does it.
version_file = 'src/toil_vcfcomp/version.py'
try:
    git_commit = subprocess.check_output(['git', 'rev-parse', 'HEAD'], encoding='ascii',
                                         stderr=subprocess.DEVNULL).strip()
except (subprocess.CalledProcessError, FileNotFoundError):
    # not a git checkout (eg installing from a tarball, or before git init)
    git_commit = 'unknown'
with open(version_file, 'w') as version_fh:
    version_fh.write("toil_vcfcomp_commit = '{}'\n".format(git_commit))

setup(
    name = 'toil-vcfcomp',
    version = '0.1.0',
    author = 'Glenn Hickey',
    description = 'Benchmark a multi-sample VCF against per-sample truth sets, in parallel with Toil',
    package_dir = {'': 'src'},
    packages = find_packages(where='src'),
    zip_safe = False,

    python_requires = '>=3.10',

    # note: toil is deliberately not here.  it lives in toil-requirement.txt so it can be
    # exact-pinned and installed first, as in cactus.
    install_requires = [
        'pysam'],

    entry_points = {
        'console_scripts': ['vcfcomp = toil_vcfcomp.vcfcomp:main',
                            'hprc-truthset = toil_vcfcomp.hprc_truthset:main',
                            'giab-truthset = toil_vcfcomp.giab_truthset:main',
                            'toil-vcfcomp = toil_vcfcomp.toil_vcfcomp:main']},)
