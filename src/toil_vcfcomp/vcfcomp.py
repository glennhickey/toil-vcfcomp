#!/usr/bin/env python3
"""
VCF comparison with hap.py and truvari, alongiside some postproecssing functions

Requirements:  * python3 environment with pysam installed
               * bcftools, samtools and tabix
               * docker (or singularity, or the tools themselves: see --binariesMode)

Example:

vcfcomp.py download-q100

vcfcomp.py happy --truth CHM13v2.0_HG2-T2TQ100-V1.1_smvar.vcf.gz \
 --calls hprc-v2.0-mc-chm13.wave.vcf.gz --ref hs1.fa \
 --regions CHM13v2.0_HG2-T2TQ100-V1.1_smvar.benchmark.bed \
 --out-dir happy-hprc-v2.0-mc-chm13.HG002.wave --sample HG002 \
 --threads 8 --exclude-y --max-length 50

vcfcomp.py happy-breakdown --vcf happy-hprc-v2.0-mc-chm13.HG002.wave/happy.output.vcf.gz




"""

import os, sys
import subprocess
import pysam
import argparse
import shlex
import shutil
from collections import defaultdict

from toil_vcfcomp.common import run_tool, add_binaries_options, setup_binaries, binaries_mode, MNT

default_happy_docker = 'jmcdani20/hap.py:v0.3.12'
default_happy_options = '--gender=male --pass-only --engine=vcfeval'
default_happy_max_length = 1000
default_truvari_docker = 'solyris/truvari:v5.3'
#these ones are from hprc v1 paper (but using --pick multi instead of --multimatch)
default_truvari_options = '-O 0.0 -r 1000 -p 0.0 -P 0.3 -C 1000 -s 50 -S 15 --sizemax 100000 --no-ref c --pick multi'
# these are from https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/analysis/NIST_HG002_DraftBenchmark_defrabbV0.019-20241113/NIST_HG002_DraftBenchmark_defrabbV0.019-20241113_README.md
#default_truvari_options = '--pick ac --passonly -r 2000 -C 5000'
default_vcfeval_docker = 'kockan/vcfeval_docker:v1.1'
default_vcfeval_options = '--decompose --ref-overlap'
default_aardvark_max_length = 100000
default_aardvark_min_length = 0
default_aardvark_docker = 'quay.io/glennhickey/aardvark:v0.9.0'
default_aardvark_options = '--enable-haplotype-metrics --enable-weighted-haplotype-metrics --enable-record-basepair-metrics'
# dipcall flags its unreliable calls with HET1/HET2/GAP1/GAP2/DIPX but leaves everything else
# FILTER="." rather than PASS.  so the filter has to keep "." too -- a --pass-only style filter
# would silently throw away about 96% of a dipcall truth set.
default_filter_expr = 'FILTER="." || FILTER="PASS"'
# GRCh38 pseudoautosomal regions.  a male is genuinely diploid here, so chrX in a PAR must not be
# flattened to haploid.  1-based inclusive, matching what dipcall uses for its -x par bed.
default_par_regions = {'chrx': [(10001, 2781479), (155701383, 156030895)]}

def run_pipe(command, stdout=None):
    """ run a shell pipeline and fail if ANY stage fails.

    A bare `sh -c 'a | b'` reports only b's exit status.  When an upstream bcftools dies partway
    through -- a truncated input, a reference missing a contig, an OOM kill -- the downstream
    stage still sees a valid header plus a prefix of records, exits 0, and leaves a perfectly
    well-formed VCF containing only part of the data.  Nothing downstream can tell the difference,
    so the sample lands in the report as an ordinary row computed against a partial truth set.
    """
    subprocess.check_call(['bash', '-o', 'pipefail', '-c', command], stdout=stdout)

def remove_vcf(vcf_path):
    """ drop an intermediate and its index together, so the .tbi does not outlive it """
    for path in [vcf_path, vcf_path + '.tbi']:
        try:
            os.remove(path)
        except OSError:
            pass

def vcf_contigs(vcf_path):
    """ set of contigs that actually carry records, from the index """
    result = subprocess.run(['bcftools', 'index', '-s', vcf_path], capture_output=True, text=True, check=True)
    return set(line.split('\t')[0] for line in result.stdout.strip().split('\n') if line.strip())

def check_contigs_preserved(input_vcf, output_vcf):
    """ post-condition for the normalize step.

    pipefail catches a stage that dies, but bcftools exits 0 on a truncated bgzf *input*, so a
    prefix can still get through cleanly.  Comparing contig sets is a cheap, index-only check that
    catches the shape of loss that matters (a whole chromosome silently disappearing).
    """
    lost = vcf_contigs(input_vcf) - vcf_contigs(output_vcf)
    if lost:
        raise RuntimeError(
            'normalization dropped every record on {} while writing {}. This usually means the '
            '--ref fasta does not contain that contig, or the input was truncated. Refusing to '
            'continue with a partial VCF.'.format(', '.join(sorted(lost)), os.path.basename(output_vcf)))

def vcf_preprocess(input_vcf,
                   output_vcf,
                   ref_fasta,
                   remove_y,
                   haploid_x,
                   sample,
                   min_length,
                   max_length,
                   bi_allelic,
                   normalize,
                   filter_expr=None,
                   par_regions=None):
    """
    hap.py will fail if chrY is diploid (ie .|1).  So we need to flatten it to haploid.  Will also extract
    A sample if given, or just remove Y entirely if specified.

    haploid_x flattens chrX outside the pseudoautosomal regions only -- inside them a male really
    does have two copies, and collapsing those would throw away every heterozygous PAR call.
    """

    output_dir = os.path.dirname(output_vcf)
    try:
        os.makedirs(output_dir)
    except:
        pass
    assert os.path.isdir(output_dir)

    tmp_vcf = output_vcf.replace('.vcf', '.c1.vcf'.format('.' + sample if sample else ''))
    view_cmd = ['bcftools', 'view', input_vcf, '-c', '1', '-Oz']
    assert max_length is not None and min_length is not None
    include_expr = 'STRLEN(REF) <= {} && STRLEN(ALT) <= {} && (STRLEN(REF) >= {} || STRLEN(ALT) >= {})'.format(max_length, max_length, min_length, min_length)
    if filter_expr:
        include_expr += ' && ({})'.format(filter_expr)
    # shlex.quote rather than wrapping in double quotes by hand: bcftools string literals are
    # themselves double quoted, so a filter like FILTER="PASS" would otherwise close the shell
    # quoting early and bcftools would see a mangled expression
    view_cmd += ['-i', shlex.quote(include_expr)]
    if sample:
        view_cmd += ['-as', sample]
    # remove tags for aardvark.  -Oz on the annotate too: without it the last stage writes plain
    # VCF, so the .vcf.gz this lands in was never actually bgzf and could not be indexed
    view_cmd = ' '.join(view_cmd) + ' | bcftools annotate -x "INFO/TYPE,INFO/LEN,INFO/ORIGIN,INFO/INV" -Oz'
    with open(tmp_vcf, 'w') as tmp_vcf_file:
        sys.stderr.write('Running vcf preprocessing: {}\n'.format(view_cmd))
        run_pipe(view_cmd, stdout=tmp_vcf_file)
    # index every intermediate, so check_contigs_preserved can read contig sets straight from the
    # index rather than rescanning the records
    subprocess.check_call(['tabix', '-fp', 'vcf', tmp_vcf])
    if bi_allelic:
        sys.stderr.write('Splitting into bi-allelic variants\n')
        norm_vcf = tmp_vcf.replace('.vcf', '.biallelic.vcf')
        with open(norm_vcf, 'w') as norm_vcf_file:
            subprocess.check_call(['bcftools', 'norm', '-m', '-any', tmp_vcf, '-Oz'],
                                  stdout=norm_vcf_file)
        subprocess.check_call(['tabix', '-fp', 'vcf', norm_vcf])
        remove_vcf(tmp_vcf)
        tmp_vcf = norm_vcf
    if normalize:
        sys.stderr.write('Normalizing indels with bcftools\n')
        norm_vcf = tmp_vcf.replace('.vcf', '.norm.vcf')
        # -T keeps bcftools sort's spill inside the output directory.  the default is $TMPDIR,
        # which under a batch system is node-local /tmp: outside the job's disk request, invisible
        # to --workDir, and left behind whenever the job is killed
        sort_tmp = os.path.join(output_dir, 'bcftools-sort-tmp')
        with open(norm_vcf, 'w') as norm_vcf_file:
            run_pipe(f'bcftools norm -f {ref_fasta} {tmp_vcf} -Oz -d exact -cw | '
                     f'bcftools sort -T {shlex.quote(sort_tmp)} | bgzip',
                     stdout=norm_vcf_file)
        shutil.rmtree(sort_tmp, ignore_errors=True)
        subprocess.check_call(['tabix', '-fp', 'vcf', norm_vcf])
        check_contigs_preserved(tmp_vcf, norm_vcf)
        remove_vcf(tmp_vcf)
        tmp_vcf = norm_vcf

    subprocess.check_call(['tabix', '-fp', 'vcf', tmp_vcf])


    tmp_vcf_file = pysam.VariantFile(tmp_vcf, 'rb')
    output_vcf_file = pysam.VariantFile(output_vcf, 'w', header=tmp_vcf_file.header)

    sys.stderr.write('Correcting/removing chrX,chrY,chrM GTs\n')
    for var in tmp_vcf_file.fetch():
        contig = var.contig.lower()
        skip = ('chry' in contig and remove_y) or 'chrm' in contig
        flatten = not skip and ('chry' in contig or
                                (haploid_x and 'chrx' in contig and not in_par(var.contig, var.pos, par_regions)))
        if flatten:
            for sample in var.samples.values():
                gt = sample['GT']
                assert len(gt) in [1,2]
                if len(gt) == 2:
                    sample['GT'] = tuple([haploid_allele(gt)])

        if not skip:
            output_vcf_file.write(var)

    tmp_vcf_file.close()
    output_vcf_file.close()
    
    subprocess.check_call(['tabix', '-fp', 'vcf', output_vcf])
    remove_vcf(tmp_vcf)

def parse_par_bed(bed_path):
    """ read a PAR bed into the {contig: [(start, end)]} shape in_par() wants.  BED is 0-based
    half-open, our intervals are 1-based inclusive, so the start shifts by one """
    regions = defaultdict(list)
    with open(bed_path) as bed_file:
        for line in bed_file:
            fields = line.split()
            if len(fields) < 3 or line.startswith(('#', 'track', 'browser')):
                continue
            regions[fields[0].lower()].append((int(fields[1]) + 1, int(fields[2])))
    if not regions:
        raise SystemExit('{} contains no usable BED intervals'.format(bed_path))
    return dict(regions)

def in_par(contig, pos, par_regions=None):
    """ is this position in a pseudoautosomal region, where a male is diploid after all """
    if par_regions is None:
        par_regions = default_par_regions
    for name, intervals in par_regions.items():
        if name in contig.lower():
            return any(start <= pos <= end for start, end in intervals)
    return False

def haploid_allele(gt):
    """ collapse a diploid genotype to the one allele that survives on a haploid contig.

    Taking gt[0] and only falling back to gt[1] when gt[0] is missing is order-dependent: it
    rescues .|1 correctly but turns the biologically identical 0|1 into a reference call while
    1|0 survives.  Prefer whichever allele is non-reference, so the collapse depends on the
    genotype rather than on which haplotype the assembler happened to put first.
    """
    called = [a for a in gt if a not in ['.', None]]
    if not called:
        return gt[0]
    non_ref = [a for a in called if a != 0]
    return non_ref[0] if non_ref else called[0]

def stage_inputs(files, output_dir):
    """ get the inputs into the output directory, since that's the only thing mounted into the
    container.  compare absolute paths -- a relative-vs-absolute mismatch used to make this try
    to link a file onto itself.  fall back to a copy when the hard link would cross a filesystem,
    which happens as soon as inputs and the job's temp dir live on different mounts """
    for f in files:
        if not os.path.isfile(f):
            continue
        if os.path.realpath(os.path.dirname(f)) == os.path.realpath(output_dir):
            continue
        dest = os.path.join(os.path.abspath(output_dir), os.path.basename(f))
        if os.path.exists(dest) or os.path.islink(dest):
            # compare inodes, not path strings: when a symlinked or bind-mounted component makes
            # the destination another name for the source, removing it destroys the input.  islink
            # is checked separately so a dangling symlink is replaced rather than written through
            if not os.path.islink(dest) and os.path.samefile(f, dest):
                continue
            os.remove(dest)
        try:
            os.link(os.path.abspath(f), dest)
        except OSError:
            shutil.copyfile(os.path.abspath(f), dest)

def stage_sdf(sdf, output_dir):
    """ get an rtg SDF into output_dir, since a directory cannot be hard-linked in.  returns the
    staged path. """
    sdf = sdf.rstrip('/')
    dest = os.path.join(os.path.abspath(output_dir), os.path.basename(sdf))
    # realpath, not abspath: if the destination is another name for the source through a symlink
    # or bind mount, the rmtree below would delete the SDF we are about to use
    if os.path.realpath(sdf) == os.path.realpath(dest):
        return dest
    if os.path.islink(dest):
        os.remove(dest)
    elif os.path.exists(dest):
        shutil.rmtree(dest)
    shutil.copytree(sdf, dest)
    return dest

def happy(truth_vcf,
          calls_vcf,
          ref_fasta,
          bed_regions,
          sample,
          output_dir,
          threads,
          options,
          docker_image,
          sdf=None):
    """
    Run hap.py on two vcfs, using the provided reference and bed.  `sdf` is a prebuilt rtg SDF:
    hap.py's default --engine=vcfeval otherwise runs rtg format over the whole reference on every
    single invocation, into scratch space nothing has accounted for.
    """

    try:
        os.makedirs(output_dir)
    except:
        pass
    assert os.path.isdir(output_dir)

    # put everything in the same place
    stage_inputs([truth_vcf, truth_vcf + '.tbi', calls_vcf, calls_vcf + '.tbi',
                  ref_fasta, ref_fasta + '.fai', bed_regions], output_dir)

    # hap.py is python2 only, so it normally comes from a container.  in local mode the binary has
    # to be found on the PATH instead of at its in-container path
    happy_bin = 'hap.py' if binaries_mode() == 'local' else '/opt/hap.py/bin/hap.py'
    happy_cmd = [happy_bin,
                 MNT + os.path.basename(truth_vcf),
                 MNT + os.path.basename(calls_vcf),
                 '-f', MNT + os.path.basename(bed_regions),
                 '-r', MNT + os.path.basename(ref_fasta),
                 '--threads', str(threads),
                 '--scratch-prefix', MNT.rstrip('/'),
                 '-o', MNT + os.path.basename(output_dir)]
    if sdf:
        happy_cmd += ['--engine-vcfeval-template', MNT + os.path.basename(stage_sdf(sdf, output_dir))]
    if options:
        happy_cmd += options.split()

    run_tool(happy_cmd, image=docker_image, work_dir=output_dir)

    # use Parsa's convention of having a happy.output.vcf.gz
    subprocess.check_call(['ln', '-sf', os.path.abspath(os.path.join(output_dir, os.path.basename(output_dir) + '.vcf.gz')),
                           os.path.abspath(os.path.join(output_dir, 'happy-output.vcf.gz'))])
    subprocess.check_call(['ln', '-sf', os.path.abspath(os.path.join(output_dir, os.path.basename(output_dir) + '.vcf.gz.tbi')),
                           os.path.abspath(os.path.join(output_dir, 'happy-output.vcf.gz.tbi'))])    

def happy_chromsplit(happy_vcf):
    """ breakdown the happy vcf into a per-chromosome table"""

    snp_fp = defaultdict(int)
    snp_fn = defaultdict(int)
    indel_fp = defaultdict(int)
    indel_fn = defaultdict(int)

    vcf_file = pysam.VariantFile(happy_vcf, 'rb')
    contigs = set(['_total_'])
    for var in vcf_file.fetch():
        assert len(var.samples.values()) == 2
        contigs.add(var.contig)
        for sample in var.samples.values():
            bd = sample['BD']
            bvt = sample['BVT']
            if bd == 'FP' and bvt == 'SNP':
                snp_fp[var.contig] += 1
                snp_fp['_total_'] += 1
            elif bd == 'FP' and bvt == 'INDEL':
                indel_fp[var.contig] += 1
                indel_fp['_total_'] += 1
            elif bd == 'FN' and bvt == 'SNP':
                snp_fn[var.contig] += 1
                snp_fn['_total_'] += 1                
            elif bd == 'FN' and bvt == 'INDEL':
                indel_fn[var.contig] += 1
                indel_fn['_total_'] += 1                
    vcf_file.close()

    output_table = [['type'] + sorted(contigs)]
    output_table += [['ERRORS'], ['SNP-FP'], ['SNP-FN'], ['INDEL-FP'], ['INDEL-FN']]
    for contig in sorted(contigs):
        output_table[-5].append(str(snp_fp[contig] + snp_fn[contig] + indel_fp[contig] + indel_fn[contig]))
        output_table[-4].append(str(snp_fp[contig]))
        output_table[-3].append(str(snp_fn[contig]))
        output_table[-2].append(str(indel_fp[contig]))
        output_table[-1].append(str(indel_fn[contig]))
    
    return output_table

def truvari(truth_vcf,
            calls_vcf,
            ref_fasta,
            bed_regions,
            sample,
            output_dir,
            options,
            docker_image):
    """
    Run truvari on two vcfs, using the provided reference and bed
    """
    try:
        os.makedirs(output_dir)
    except:
        pass
    assert os.path.isdir(output_dir)

    # truvari needs a non-existent output directory
    tv_dir = os.path.join(output_dir, 'tv')
    if os.path.exists(tv_dir):
        shutil.rmtree(tv_dir)

    # put everything in the same place
    stage_inputs([truth_vcf, truth_vcf + '.tbi', calls_vcf, calls_vcf + '.tbi',
                  ref_fasta, ref_fasta + '.fai', bed_regions], output_dir)

    # run truvari with docker (to be consistent with hap.py method, tho truvari easy to install with pip)
    truvari_cmd = ['truvari', 'bench',
                   '-c', MNT + os.path.basename(calls_vcf),
                   '-b', MNT + os.path.basename(truth_vcf),
                   '--includebed', MNT + os.path.basename(bed_regions),
                   '-f', MNT + os.path.basename(ref_fasta),
                   '-o', MNT + 'tv/']
    if options:
        truvari_cmd += options.split()

    run_tool(truvari_cmd, image=docker_image, work_dir=output_dir)

    subprocess.check_call('mv {}/* {}'.format(tv_dir, output_dir), shell=True)
    shutil.rmtree(tv_dir)

def truvari_chromsplit(truvari_outdir):
    """ breakdown the truvari vcfs into a per-chromosome table"""

    fp = defaultdict(int)
    fn = defaultdict(int)

    fp_vcf = os.path.join(truvari_outdir, 'fp.vcf.gz')
    fp_vcf_file = pysam.VariantFile(fp_vcf, 'r')
    contigs = set(['_total_'])
    for var in fp_vcf_file.fetch():
        contigs.add(var.contig)
        fp[var.contig] += 1
        fp['_total_'] += 1
    fp_vcf_file.close()

    fn_vcf = os.path.join(truvari_outdir, 'fn.vcf.gz')
    fn_vcf_file = pysam.VariantFile(fn_vcf, 'r')
    # do not reset contigs here: it has to accumulate over both passes, otherwise any contig
    # that only appears in the FP vcf is dropped from the table's columns
    for var in fn_vcf_file.fetch():
        contigs.add(var.contig)
        fn[var.contig] += 1
        fn['_total_'] += 1
    fn_vcf_file.close()

    output_table = [['type'] + sorted(contigs)]
    output_table += [['ERRORS'], ['SV-FP'], ['SV-FN']]
    for contig in sorted(contigs):
        output_table[-3].append(str(fp[contig] + fn[contig]))
        output_table[-2].append(str(fp[contig]))
        output_table[-1].append(str(fn[contig]))
    
    return output_table

def vcfeval(truth_vcf,
            calls_vcf,
            ref_fasta,
            bed_regions,
            sample,
            output_dir,
            threads,
            options,
            docker_image,
            sdf=None):
    """
    Run rtg vcfeval on two vcfs, using the provided reference and bed.  `sdf` is an already-built
    rtg SDF directory: pass one to skip the rtg format step, which is otherwise repeated on every
    single invocation and dominates the runtime when comparing many samples.
    """

    try:
        os.makedirs(output_dir)
    except:
        pass
    assert os.path.isdir(output_dir)

    # vcfeval needs a non-existent output directory
    ve_dir = os.path.join(output_dir, 've')
    if os.path.exists(ve_dir):
        shutil.rmtree(ve_dir)

    # put everything in the same place
    stage_inputs([truth_vcf, truth_vcf + '.tbi', calls_vcf, calls_vcf + '.tbi',
                  ref_fasta, ref_fasta + '.fai', bed_regions], output_dir)

    if sdf is None and ref_fasta.endswith('.SDF'):
        sdf = ref_fasta
    if sdf:
        ref_name = os.path.basename(stage_sdf(sdf, output_dir))
    else:
        ref_name = 'ref.SDF'
        if os.path.exists(os.path.join(output_dir, ref_name)):
            shutil.rmtree(os.path.join(output_dir, ref_name))
        run_tool(['rtg', 'format', MNT + os.path.basename(ref_fasta), '-o', MNT + ref_name],
                 image=docker_image, work_dir=output_dir)

    # run vcfeval
    vcfeval_cmd = ['rtg', 'vcfeval',
                   '-b', MNT + os.path.basename(truth_vcf),
                   '-c', MNT + os.path.basename(calls_vcf),                 
                   '-e', MNT + os.path.basename(bed_regions),
                   '-t', MNT + ref_name,
                   '--threads', str(threads),
                   '-o', MNT + os.path.basename(ve_dir)]
    if sample:
        vcfeval_cmd += ['--sample', sample]
    if options:
        vcfeval_cmd += options.split()
    run_tool(vcfeval_cmd, image=docker_image, work_dir=output_dir)

    subprocess.check_call('mv {}/* {}'.format(ve_dir, output_dir), shell=True)
    shutil.rmtree(ve_dir)

def vcfeval_chromsplit(truvari_outdir):
    """ breakdown the truvari vcfs into a per-chromosome table"""

    snp_fp = defaultdict(int)
    snp_fn = defaultdict(int)
    indel_fp = defaultdict(int)
    indel_fn = defaultdict(int)

    def is_snp(var):
        """ note: this isn't very sophisticated and may differe from hap.py's sense quite a bit """
        return all(len(a) == len(var.alleles[0]) for a in var.alleles)
    
    fp_vcf = os.path.join(truvari_outdir, 'fp.vcf.gz')
    fp_vcf_file = pysam.VariantFile(fp_vcf, 'r')
    contigs = set(['_total_'])
    for var in fp_vcf_file.fetch():        
        contigs.add(var.contig)
        if is_snp(var):
            snp_fp[var.contig] += 1
            snp_fp['_total_'] += 1
        else:
            indel_fp[var.contig] += 1
            indel_fp['_total_'] += 1
    fp_vcf_file.close()

    fn_vcf = os.path.join(truvari_outdir, 'fn.vcf.gz')
    fn_vcf_file = pysam.VariantFile(fn_vcf, 'r')
    # do not reset contigs here: it has to accumulate over both passes, otherwise any contig
    # that only appears in the FP vcf is dropped from the table's columns
    for var in fn_vcf_file.fetch():
        contigs.add(var.contig)
        if is_snp(var):
            snp_fn[var.contig] += 1
            snp_fn['_total_'] += 1
        else:
            indel_fn[var.contig] += 1
            indel_fn['_total_'] += 1
    fn_vcf_file.close()

    output_table = [['type'] + sorted(contigs)]
    output_table += [['ERRORS'], ['SNP-FP'], ['SNP-FN'], ['INDEL-FP'], ['INDEL-FN']]
    for contig in sorted(contigs):
        output_table[-5].append(str(snp_fp[contig] + snp_fn[contig] + indel_fp[contig] + indel_fn[contig]))
        output_table[-4].append(str(snp_fp[contig]))
        output_table[-3].append(str(snp_fn[contig]))
        output_table[-2].append(str(indel_fp[contig]))
        output_table[-1].append(str(indel_fn[contig]))
    
    return output_table

def aardvark(truth_vcf,
             calls_vcf,
             ref_fasta,
             bed_regions,
             sample,
             output_dir,
             threads,
             options,
             docker_image):
    """
    Run Aardvark to do the vcf comparison
    """
    try:
        os.makedirs(output_dir)
    except:
        pass
    assert os.path.isdir(output_dir)

    # fresh output directory
    ve_dir = os.path.join(output_dir, 'ard')
    if os.path.exists(ve_dir):
        shutil.rmtree(ve_dir)

    # put everything in the same place
    stage_inputs([truth_vcf, truth_vcf + '.tbi', calls_vcf, calls_vcf + '.tbi',
                  ref_fasta, ref_fasta + '.fai', bed_regions], output_dir)

    # run aardvark
    aardvark_cmd = ['aardvark', 'compare',
                    '-t', MNT + os.path.basename(truth_vcf),
                    '-q', MNT + os.path.basename(calls_vcf),                 
                    '-b', MNT + os.path.basename(bed_regions),
                    '-r', MNT + os.path.basename(ref_fasta),
                    '--threads', str(threads),
                    '-o', MNT + os.path.basename(ve_dir),
                    '--compare-label', os.path.basename(calls_vcf)]
    if options:
        aardvark_cmd += options.split()
    run_tool(aardvark_cmd, image=docker_image, work_dir=output_dir)

    subprocess.check_call('mv {}/* {}'.format(ve_dir, output_dir), shell=True)
    shutil.rmtree(ve_dir)
    
def aardvark_chromsplit(aardvark_outdir):
    """ breakdown the aardvark vcfs into a per-chromosome table"""

    snp_fp = defaultdict(int)
    snp_fn = defaultdict(int)
    indel_fp = defaultdict(int)
    indel_fn = defaultdict(int)
    sv_snp_fp = defaultdict(int)
    sv_snp_fn = defaultdict(int)
    sv_indel_fp = defaultdict(int)
    sv_indel_fn = defaultdict(int)
    
    def is_snp(var):
        """ note: this isn't very sophisticated and may differe from hap.py's sense quite a bit """
        return all(len(a) == len(var.alleles[0]) for a in var.alleles)
    def is_sv(var):
        """ this either """
        return any(len(a) >= 50 for a in var.alleles)
    
    fp_vcf = os.path.join(aardvark_outdir, 'query.vcf.gz')
    fp_vcf_file = pysam.VariantFile(fp_vcf, 'r')
    contigs = set(['_total_'])
    for var in fp_vcf_file.fetch():
        assert len(var.samples.values()) == 1
        sample = var.samples.values()[0]
        if sample['BD'] == 'FP':            
            contigs.add(var.contig)
            if is_sv(var):
                snp_dict = sv_snp_fp
                indel_dict = sv_indel_fp
            else:
                snp_dict = snp_fp
                indel_dict = indel_fp
            if is_snp(var):
                snp_dict[var.contig] += 1
                snp_dict['_total_'] += 1
            else:
                indel_dict[var.contig] += 1
                indel_dict['_total_'] += 1
    fp_vcf_file.close()

    fn_vcf = os.path.join(aardvark_outdir, 'truth.vcf.gz')
    fn_vcf_file = pysam.VariantFile(fn_vcf, 'r')
    # do not reset contigs here: it has to accumulate over both passes, otherwise any contig
    # that only appears in the FP vcf is dropped from the table's columns
    for var in fn_vcf_file.fetch():
        assert len(var.samples.values()) == 1
        sample = var.samples.values()[0]
        if sample['BD'] == 'FN':            
            contigs.add(var.contig)
            if is_sv(var):
                snp_dict = sv_snp_fn
                indel_dict = sv_indel_fn
            else:
                snp_dict = snp_fn
                indel_dict = indel_fn
            if is_snp(var):
                snp_dict[var.contig] += 1
                snp_dict['_total_'] += 1
            else:
                indel_dict[var.contig] += 1
                indel_dict['_total_'] += 1
    fn_vcf_file.close()

    output_table = [['type'] + sorted(contigs)]
    output_table += [['ERRORS'], ['SNP-FP'], ['SNP-FN'], ['INDEL-FP'], ['INDEL-FN'],
                     ['SV-SNP-FP'], ['SV-SNP-FN'], ['SV-INDEL-FP'], ['SV-INDEL-FN']]
    for contig in sorted(contigs):
        output_table[-9].append(str(snp_fp[contig] + snp_fn[contig] + indel_fp[contig] + indel_fn[contig] +
                                    sv_snp_fp[contig] + sv_snp_fn[contig] + sv_indel_fp[contig] + sv_indel_fn[contig]))
        output_table[-8].append(str(snp_fp[contig]))
        output_table[-7].append(str(snp_fn[contig]))
        output_table[-6].append(str(indel_fp[contig]))
        output_table[-5].append(str(indel_fn[contig]))
        output_table[-4].append(str(sv_snp_fp[contig]))
        output_table[-3].append(str(sv_snp_fn[contig]))
        output_table[-2].append(str(sv_indel_fp[contig]))
        output_table[-1].append(str(sv_indel_fn[contig]))        
    
    return output_table

def download_q100(dict_only=False):
    """
    Get all the benchmark data needed for the draft HG002 Q100 benchmarks
    """
    base_url = 'https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/analysis/NIST_HG002_DraftBenchmark_defrabbV0.019-20241113/'

    name_dict = {}
    for ref in ['CHM13v2.0', 'GRCh38']:
        for vartype in ['smvar', 'stvar']:
            for ext in ['benchmark.bed', 'vcf.gz', 'vcf.gz.tbi']:
                name = '{}_HG2-T2TQ100-V1.1_{}.{}'.format(ref, vartype, ext)
                if not dict_only:
                    sys.stderr.write('Downloading {}\n'.format(name))
                    subprocess.check_call(['wget', '-q', os.path.join(base_url, name), '-O', name])
                    # hack to drop y from chm13 benchmark
                    if ref == 'CHM13v2.0' and ext == 'benchmark.bed':
                        no_y_name = name.replace('.benchmark.bed', '.no.y.benchmark.bed')
                        with open(no_y_name, 'w') as no_y_bed_file:
                            subprocess.check_call(['grep', '-v', '^chrY', name], stdout=no_y_bed_file)
                        name = no_y_name
                name_dict[(ref, vartype, ext)] = name

    if not dict_only:
        sys.stderr.write('Downlading hs1.fa.gz\n')
        subprocess.check_call(['wget', '-q', 'https://hgdownload.soe.ucsc.edu/goldenPath/hs1/bigZips/hs1.fa.gz', '-O', 'hs1.fa.gz'])

        sys.stderr.write('Downloading hg38.fa.gz\n')
        subprocess.check_call(['wget', '-q', 'https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz', '-O', 'hg38.fa.gz'])

    for ref in 'hs1', 'hg38':
        # note: I'm pretty sure hap.py doesn't work properly with compressed references
        if not dict_only:
            sys.stderr.write('Indexing {}.fa\n'.format(ref))
            subprocess.check_call(['gzip', '-fd', '{}.fa.gz'.format(ref)])
            subprocess.check_call(['samtools', 'faidx', '{}.fa'.format(ref)])
        name_dict[ref] = '{}.fa'.format(ref)

        #hack
        #name_dict['hg38'] = 'GCA_000001405.15_GRCh38_no_alt_analysis_set_maskedGRC_exclusions_v2.fasta'
    return name_dict
                                 
def eval_q100(out_dir, grch38_vcfs, chm13_vcfs, download, min_length, max_length, threads, eval_type,
              happy_options, truvari_options, vcfeval_options, aardvark_options):
    """ run q100 hap.py and truvari evaluation on set of input vcfs and
    tabulate the results """

    name_dict = download_q100(dict_only = not download)

    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    assert os.path.isdir(out_dir)

    results_table = {}

    for i, vcf in enumerate(grch38_vcfs + chm13_vcfs):
        is_grch38 = i < len(grch38_vcfs)
        happy_out_dir = os.path.join(out_dir, 'happy-' + os.path.basename(vcf.replace('.vcf.gz','')))
        # run hap.py preprocessing
        hap_calls = os.path.join(happy_out_dir, os.path.basename(vcf.replace('.vcf.gz', '.hap.vcf.gz')))
        if 'happy' in eval_type:
            vcf_preprocess(vcf,
                           hap_calls,
                           name_dict['hg38' if is_grch38 else 'hs1'],
                           not is_grch38,
                           True,
                           'HG002',
                           min_length,
                           max_length,
                           False,
                           True)        
            # run hap.py
            happy(name_dict[('GRCh38' if is_grch38 else 'CHM13v2.0', 'smvar', 'vcf.gz')],
                  hap_calls,
                  name_dict['hg38' if is_grch38 else 'hs1'],
                  name_dict[('GRCh38' if is_grch38 else 'CHM13v2.0', 'smvar', 'benchmark.bed')],
                  'HG002',
                  happy_out_dir,
                  threads,
                  happy_options,
                  default_happy_docker)

            # run hap.py postprocessing
            results_table[os.path.basename(happy_out_dir)] = happy_chromsplit(os.path.join(happy_out_dir, 'happy-output.vcf.gz'))   

        if 'truvari' in eval_type:
            truvari_out_dir = os.path.join(out_dir, 'truvari-' + os.path.basename(vcf.replace('.vcf.gz','')))
            # run truvari preprocessing
            tru_calls = os.path.join(truvari_out_dir, os.path.basename(vcf.replace('.vcf.gz', '.tru.vcf.gz')))
            vcf_preprocess(vcf,
                           tru_calls,
                           name_dict['hg38' if is_grch38 else 'hs1'],
                           not is_grch38,
                           True,
                           'HG002',
                           min_length,
                           max_length,
                           True,
                           True)
        
            # run truvari
            truvari(name_dict[('GRCh38' if is_grch38 else 'CHM13v2.0', 'stvar', 'vcf.gz')],
                    tru_calls,
                    name_dict['hg38' if is_grch38 else 'hs1'],
                    name_dict[('GRCh38' if is_grch38 else 'CHM13v2.0', 'stvar', 'benchmark.bed')],
                    'HG002',
                    truvari_out_dir,
                    truvari_options,
                    default_truvari_docker)

            # run truvari postprocessing
            results_table[os.path.basename(truvari_out_dir)] = truvari_chromsplit(truvari_out_dir)

        if 'vcfeval' in eval_type:
            vcfeval_out_dir = os.path.join(out_dir, 'vcfeval-' + os.path.basename(vcf.replace('.vcf.gz','')))
            # run vcfeval preprocessing
            hap_calls = os.path.join(vcfeval_out_dir, os.path.basename(vcf.replace('.vcf.gz', '.hap.vcf.gz')))            
            vcf_preprocess(vcf,
                           hap_calls,
                           name_dict['hg38' if is_grch38 else 'hs1'],
                           not is_grch38,
                           True,
                           'HG002',
                           min_length,
                           max_length,
                           False,
                           True)
            # run vcfeval
            vcfeval(name_dict[('GRCh38' if is_grch38 else 'CHM13v2.0', 'smvar', 'vcf.gz')],
                    hap_calls,
                    name_dict['hg38' if is_grch38 else 'hs1'],
                    name_dict[('GRCh38' if is_grch38 else 'CHM13v2.0', 'smvar', 'benchmark.bed')],
                    'HG002',
                    vcfeval_out_dir,
                    threads,
                    vcfeval_options,
                    default_vcfeval_docker)

            # run vcfeval postprocessing
            results_table[os.path.basename(vcfeval_out_dir)] = vcfeval_chromsplit(vcfeval_out_dir)

        if 'aardvark' in eval_type:
            aardvark_out_dir = os.path.join(out_dir, 'aardvark-' + os.path.basename(vcf.replace('.vcf.gz','')))
            # run aardvark preprocessing
            hap_calls = os.path.join(aardvark_out_dir, os.path.basename(vcf.replace('.vcf.gz', '.hap.vcf.gz')))
            vcf_preprocess(vcf,
                           hap_calls,
                           name_dict['hg38' if is_grch38 else 'hs1'],
                           not is_grch38,
                           True,
                           'HG002',
                           min_length,
                           default_aardvark_max_length,
                           False,
                           True)
            # run aardvark
            # note we compare all variants here
            aardvark(name_dict[('GRCh38' if is_grch38 else 'CHM13v2.0', 'stvar', 'vcf.gz')],
                     hap_calls,
                     name_dict['hg38' if is_grch38 else 'hs1'],
                     name_dict[('GRCh38' if is_grch38 else 'CHM13v2.0', 'stvar', 'benchmark.bed')],
                     'HG002',
                     aardvark_out_dir,
                     threads,
                     aardvark_options,
                     default_aardvark_docker)

            # run aardvark postprocessing
            results_table[os.path.basename(aardvark_out_dir)] = aardvark_chromsplit(aardvark_out_dir)

    for et in eval_type:
        with open(os.path.join(out_dir, '{}.tsv'.format(et)), 'w') as et_file:
            for k,v in results_table.items():
                if et in k:
                    for line in v:
                        et_file.write('{}\t{}\n'.format(k, '\t'.join(line)))                        
    
def main(command_line=None):                     
    parser = argparse.ArgumentParser('VCF comparison tools for evaluating pangenome graphs')
    add_binaries_options(parser)
    subparsers = parser.add_subparsers(dest='command')

    happy_parser = subparsers.add_parser('happy', help='Run hap.py')
    happy_parser.add_argument('--truth', required=True,
                       help='true VCF')
    happy_parser.add_argument('--calls', required=True,
                       help='calls VCF')
    happy_parser.add_argument('--ref', required=True,
                       help='reference fasta')
    happy_parser.add_argument('--regions', required=True,
                       help='BED regions to evaluate')
    happy_parser.add_argument('--out-dir', required=True,
                       help='output directory')    
    happy_parser.add_argument('--sample',
                       help='subset to this sample')
    happy_parser.add_argument('--max-length', type=int, default=default_happy_max_length,
                       help='ignore sites with alleles longer than this')    
    happy_parser.add_argument('--docker', default=default_happy_docker,
                       help='use this docker image instead of the default')
    happy_parser.add_argument('--options', default=default_happy_options,
                       help='use these options instead of the default (surround in quotes)')
    happy_parser.add_argument('--threads', type=int, default=8,
                       help='number of threads (default=8)')
    happy_parser.add_argument('--haploid-x', action='store_true',
                              help='make sure chrX is haploid')    
    happy_parser.add_argument('--exclude-y', action='store_true',
                       help='completely ignore chrY')

    happy_breakdown_parser = subparsers.add_parser('happy-breakdown', help='Make chromosome-decomposed table of hap.py results')
    happy_breakdown_parser.add_argument('--vcf', required=True,
                                        help='happy output VCF')

    truvari_parser = subparsers.add_parser('truvari', help='Run truvari')
    truvari_parser.add_argument('--truth', required=True,
                       help='true VCF')
    truvari_parser.add_argument('--calls', required=True,
                       help='calls VCF')
    truvari_parser.add_argument('--ref', required=True,
                       help='reference fasta')
    truvari_parser.add_argument('--regions', required=True,
                       help='BED regions to evaluate')
    truvari_parser.add_argument('--out-dir', required=True,
                       help='output directory')    
    truvari_parser.add_argument('--sample',
                       help='subset to this sample')
    truvari_parser.add_argument('--docker', default=default_truvari_docker,
                                help='use this docker image instead of the default')
    truvari_parser.add_argument('--options', default=default_truvari_options,
                                help='use these options instead of the default (surround in quotes)')
    truvari_parser.add_argument('--haploid-x', action='store_true',
                              help='make sure chrX is haploid')    
    truvari_parser.add_argument('--exclude-y', action='store_true',
                              help='completely ignore chrY')

    truvari_breakdown_parser = subparsers.add_parser('truvari-breakdown', help='Make chromosome-decomposed table of truvari results')
    truvari_breakdown_parser.add_argument('--dir', required=True,
                                        help='truvari output directory')

    vcfeval_parser = subparsers.add_parser('vcfeval', help='Run rtg vcfeval')
    vcfeval_parser.add_argument('--truth', required=True,
                       help='true VCF')
    vcfeval_parser.add_argument('--calls', required=True,
                       help='calls VCF')
    vcfeval_parser.add_argument('--ref', required=True,
                       help='reference fasta (or save time by passing in a .SDF from rtg format)')
    vcfeval_parser.add_argument('--sdf',
                       help='prebuilt rtg SDF directory for --ref, to skip the rtg format step')
    vcfeval_parser.add_argument('--regions', required=True,
                       help='BED regions to evaluate')
    vcfeval_parser.add_argument('--out-dir', required=True,
                       help='output directory')
    vcfeval_parser.add_argument('--max-length', type=int, default=default_happy_max_length,
                                help='ignore sites with alleles longer than this')        
    vcfeval_parser.add_argument('--sample',
                       help='subset to this sample')
    vcfeval_parser.add_argument('--docker', default=default_vcfeval_docker,
                                help='use this docker image instead of the default')
    vcfeval_parser.add_argument('--threads', type=int, default=8,
                                help='number of threads (default=8)')    
    vcfeval_parser.add_argument('--options', default=default_vcfeval_options,
                                help='use these options instead of the default (surround in quotes)')
    vcfeval_parser.add_argument('--haploid-x', action='store_true',
                              help='make sure chrX is haploid')    
    vcfeval_parser.add_argument('--exclude-y', action='store_true',
                              help='completely ignore chrY')
    
    vcfeval_breakdown_parser = subparsers.add_parser('vcfeval-breakdown', help='Make chromosome-decomposed table of vcfeval results')
    vcfeval_breakdown_parser.add_argument('--dir', required=True,
                                        help='vcfeval output directory')

    aardvark_parser = subparsers.add_parser('aardvark', help='Run aardvark')
    aardvark_parser.add_argument('--truth', required=True,
                                 help='true VCF')
    aardvark_parser.add_argument('--calls', required=True,
                                 help='calls VCF')
    aardvark_parser.add_argument('--ref', required=True,
                                 help='reference fasta (or save time by passing in a .SDF from rtg format)')
    aardvark_parser.add_argument('--regions', required=True,
                                 help='BED regions to evaluate')
    aardvark_parser.add_argument('--out-dir', required=True,
                                 help='output directory')
    aardvark_parser.add_argument('--min-length', type=int, default=default_aardvark_min_length,
                                help='ignore sites with alleles shorter than this')
    aardvark_parser.add_argument('--max-length', type=int, default=default_aardvark_max_length,
                                help='ignore sites with alleles longer than this')            
    aardvark_parser.add_argument('--sample',
                                 help='subset to this sample')
    aardvark_parser.add_argument('--docker', default=default_aardvark_docker,
                                 help='use this docker image instead of the default')
    aardvark_parser.add_argument('--threads', type=int, default=8,
                                 help='number of threads (default=8)')    
    aardvark_parser.add_argument('--options', default=default_aardvark_options,
                                 help='use these options instead of the default (surround in quotes)')
    aardvark_parser.add_argument('--haploid-x', action='store_true',
                                 help='make sure chrX is haploid')    
    aardvark_parser.add_argument('--exclude-y', action='store_true',
                                 help='completely ignore chrY')

    aardvark_breakdown_parser = subparsers.add_parser('aardvark-breakdown', help='Make chromosome-decomposed table of aardvark results')
    aardvark_breakdown_parser.add_argument('--dir', required=True,
                                           help='aardvark output directory')
    
    hg002_q100_download_parser = subparsers.add_parser('download-q100', help='Download all HG002 t2t-q100 truth set files')

    hg002_q100_eval_parser = subparsers.add_parser('eval-q100', help='Run t2t-q100 HG002 evaluation on some graphs')
    hg002_q100_eval_parser.add_argument('--download', action='store_true',
                                        help='automatically download benchmark data into current directory')
    hg002_q100_eval_parser.add_argument('--out-dir',  required=True,
                                        help='automatically download benchmark data into current directory')
    hg002_q100_eval_parser.add_argument('--min-length', type=int, default=default_aardvark_min_length,
                                        help='ignore sites with alleles shorter than this')    
    hg002_q100_eval_parser.add_argument('--max-length', type=int, default=default_happy_max_length,
                                        help='ignore sites with alleles longer than this (does not apply to truvari)')
    hg002_q100_eval_parser.add_argument('--threads', type=int, default=8,
                                        help='number of threads (does not apply to truvari) (default=8)')    
    hg002_q100_eval_parser.add_argument('--grch38-vcfs', nargs='*', default=[],
                                        help='vcfs for HG002-T2T-Q100-GRCh38 evaluation')
    hg002_q100_eval_parser.add_argument('--chm13-vcfs', nargs='*', default=[],
                                        help='vcfs for HG002-T2T-Q100-GRCh38 evaluation')
    hg002_q100_eval_parser.add_argument('--eval-type', nargs='+', default=['vcfeval', 'truvari'],
                                        choices = ['happy', 'truvari', 'vcfeval', 'aardvark', 'aardvark-sv'],
                                        help='evaluation to run in {happy, truvari, vcfeval, aardvark, aardvark-sv}, multiple allowed')
    hg002_q100_eval_parser.add_argument('--happy-options', default=default_happy_options,
                                        help='use these options instead of the default (surround in quotes)')    
    hg002_q100_eval_parser.add_argument('--truvari-options', default=default_truvari_options,
                                        help='use these options instead of the default (surround in quotes)')
    hg002_q100_eval_parser.add_argument('--vcfeval-options', default=default_vcfeval_options,
                                        help='use these options instead of the default (surround in quotes)')
    hg002_q100_eval_parser.add_argument('--aardvark-options', default=default_aardvark_options,
                                        help='use these options instead of the default (surround in quotes)')        

    # shared by the four subcommands that run a comparison
    for tool_parser in [happy_parser, truvari_parser, vcfeval_parser, aardvark_parser]:
        tool_parser.add_argument('--par-bed',
                                 help='BED of pseudoautosomal regions, where --haploid-x must not '
                                 'flatten genotypes. defaults to the GRCh38 PARs; pass CHM13 ones here')
    happy_parser.add_argument('--sdf',
                              help='prebuilt rtg SDF, passed to hap.py as --engine-vcfeval-template '
                              'so it does not rebuild one per run')

    args = parser.parse_args(command_line)
    par_regions = parse_par_bed(args.par_bed) if getattr(args, 'par_bed', None) else None

    # only the commands that actually run a tool need this: no point probing for docker just to
    # tabulate an output directory
    if args.command in ['happy', 'truvari', 'vcfeval', 'aardvark', 'eval-q100']:
        setup_binaries(args)

    if args.command in ['happy', 'truvari', 'vcfeval', 'aardvark']:
        if not os.path.isdir(args.out_dir):
            os.makedirs(args.out_dir)

        assert args.truth.endswith('.vcf.gz')
        assert args.calls.endswith('.vcf.gz')
        assert args.truth != args.calls
        assert os.path.isfile(args.truth) and os.path.isfile(args.truth + '.tbi')
        assert os.path.isfile(args.calls) and os.path.isfile(args.calls + '.tbi')
        assert os.path.isfile(args.ref) and os.path.isfile(args.regions)

    if args.command == 'happy':
        assert not args.ref.endswith('.gz')
        # run some preprocessing (only on calls -- assume truth from giab is ready to go)
        hap_calls = os.path.join(args.out_dir, os.path.basename(args.calls.replace('.vcf.gz', '.hap.vcf.gz')))        
        vcf_preprocess(args.calls,
                       hap_calls,
                       args.ref,
                       args.exclude_y,
                       args.haploid_x,
                       args.sample,
                       0,
                       args.max_length,
                       True,
                       True,
                       par_regions=par_regions)

        # run happy
        happy(args.truth, hap_calls, args.ref, args.regions, args.sample, args.out_dir,
              threads = args.threads,
              options=args.options,
              docker_image = args.docker,
              sdf = args.sdf)

    elif args.command == 'truvari':
        # run some preprocessing (only on calls -- assume truth from giab is ready to go)
        tv_calls = os.path.join(args.out_dir, os.path.basename(args.calls.replace('.vcf.gz', '.tv.vcf.gz')))        
        vcf_preprocess(args.calls,
                       tv_calls,
                       args.ref,
                       args.exclude_y,
                       args.haploid_x,
                       args.sample,
                       default_aardvark_min_length,
                       default_aardvark_max_length,
                       True,
                       True,
                       par_regions=par_regions)

        # run truvari
        truvari(args.truth, tv_calls, args.ref, args.regions, args.sample, args.out_dir,
                options=args.options,
                docker_image = args.docker)

    if args.command == 'vcfeval':
        assert not args.ref.endswith('.gz')
        # run some preprocessing (only on calls -- assume truth from giab is ready to go)
        hap_calls = os.path.join(args.out_dir, os.path.basename(args.calls.replace('.vcf.gz', '.hap.vcf.gz')))        
        vcf_preprocess(args.calls,
                       hap_calls,
                       args.ref,
                       args.exclude_y,
                       args.haploid_x,
                       args.sample,
                       0,
                       args.max_length,
                       False,
                       True,
                       par_regions=par_regions)

        # run vcfeval
        vcfeval(args.truth, hap_calls, args.ref, args.regions, args.sample, args.out_dir,
                threads = args.threads,
                options=args.options,
                docker_image = args.docker,
                sdf = args.sdf)

    if args.command == 'aardvark':
        assert not args.ref.endswith('.gz')
        # run some preprocessing (only on calls -- assume truth from giab is ready to go)
        hap_calls = os.path.join(args.out_dir, os.path.basename(args.calls.replace('.vcf.gz', '.hap.vcf.gz')))        
        vcf_preprocess(args.calls,
                       hap_calls,
                       args.ref,
                       args.exclude_y,
                       args.haploid_x,
                       args.sample,
                       args.min_length,
                       args.max_length,
                       False,
                       True,
                       par_regions=par_regions)
        
        # run aardvark
        aardvark(args.truth, hap_calls, args.ref, args.regions, args.sample, args.out_dir,
                 threads = args.threads,
                 options=args.options,
                 docker_image = args.docker)        

    elif args.command == 'happy-breakdown':
        table = happy_chromsplit(args.vcf)
        for row in table:
            print('\t'.join(row))

    elif args.command == 'truvari-breakdown':
        table = truvari_chromsplit(args.dir)
        for row in table:
            print('\t'.join(row))

    elif args.command == 'vcfeval-breakdown':
        table = vcfeval_chromsplit(args.dir)
        for row in table:
            print('\t'.join(row))

    elif args.command == 'aardvark-breakdown':
        table = aardvark_chromsplit(args.dir)
        for row in table:
            print('\t'.join(row))            
            
    elif args.command == 'download-q100':
        download_q100()

    elif args.command == 'eval-q100':
        eval_q100(args.out_dir,
                  args.grch38_vcfs,
                  args.chm13_vcfs,
                  args.download,
                  args.min_length,
                  args.max_length,
                  args.threads,
                  args.eval_type,
                  args.happy_options,
                  args.truvari_options,
                  args.vcfeval_options,
                  args.aardvark_options)

            
if __name__ == '__main__':
    main()
