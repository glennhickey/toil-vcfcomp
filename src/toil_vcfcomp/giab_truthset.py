#!/usr/bin/env python3
"""
Turn the GIAB benchmark FTP tree into a truth manifest for toil-vcfcomp.

GIAB publishes seven benchmark samples (HG001-HG007), and for HG002 a good deal more than the
others: v4.2.1 small variants, the medical-gene benchmark CMRG, the SV Tier1 set, tandem repeats,
chrXY, and the draft Q100 (v5.0q).  Each becomes a `truthset` row, so one toil-vcfcomp run can
score a call set against all of them at once and report them side by side.

The tree is NOT uniformly shaped -- some benchmarks put files directly in the benchmark directory,
some under a per-reference subdirectory, some under a further SmallVariant/StructuralVariant
split -- so this walks each benchmark directory rather than composing a path template.  Within a
directory, a VCF is paired with the BED whose name shares the longest prefix with it, which is what
actually holds across every layout GIAB uses:

  HG002_GRCh38_CMRG_smallvar_v1.00.vcf.gz   <-> HG002_GRCh38_CMRG_smallvar_v1.00.bed
  HG002_GRCh38_v5.0q_smvar.vcf.gz           <-> HG002_GRCh38_v5.0q_smvar.benchmark.bed
  HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz <-> HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed

Usage:

  giab-truthset -o giab.index.csv                          # v4.2.1 for all 7 samples, GRCh38
  giab-truthset --benchmarks all -o giab.index.csv         # every benchmark set there is
  giab-truthset --samples HG002 --benchmarks all --ref CHM13v2.0 -o giab.chm13.csv

Note that the benchmarks differ in what they cover: v4.2.1 and CMRG smallvar are small variants
(pair them with vcfeval or hap.py), SV Tier1 and CMRG SV are structural (truvari, aardvark).
Running a small-variant tool against an SV benchmark will "work" and report nonsense.

Python 3, stdlib only.
"""

import os, sys
import argparse
import re
import urllib.parse

from toil_vcfcomp.manifest import (manifest_header, write_manifest, download_rows,
                                   parse_sample_list, list_html_dir)

default_root = 'https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release'
default_ref = 'GRCh38'
default_benchmarks = 'NISTv4.2.1'

# stable, and small enough to state outright rather than scrape
giab_samples = {'HG001': 'NA12878_HG001',
                'HG002': 'AshkenazimTrio/HG002_NA24385_son',
                'HG003': 'AshkenazimTrio/HG003_NA24149_father',
                'HG004': 'AshkenazimTrio/HG004_NA24143_mother',
                'HG005': 'ChineseTrio/HG005_NA24631_son',
                'HG006': 'ChineseTrio/HG006_NA24694_father',
                'HG007': 'ChineseTrio/HG007_NA24695_mother'}

known_refs = ['GRCh38', 'GRCh37', 'CHM13v2.0', 'CHM13v1.0']

# benchmarks whose filenames carry no reference token and which live in no per-reference
# subdirectory.  there is nothing in the tree to infer these from, so they are stated.
benchmark_ref_overrides = {'NIST_SV_v0.6': 'GRCh37'}

# directories under a benchmark that hold supporting material rather than benchmarks
skip_dirs = set(['SupplementaryFiles', 'GIAB_Evaluations', 'defrabb_files', 'dipcall_output',
                 'hifiasm-assembly', 'supplementaryFiles'])

# superseded releases; --benchmarks all should not drag these in
legacy_benchmark_re = re.compile(r'^(NISTv[123]|GIABPedigree|NISTv4\.1$|latest$)')


def sample_url(root, sample):
    return '{}/{}'.format(root.rstrip('/'), giab_samples[sample])


def list_benchmarks(root, sample):
    """ benchmark set directories published for one sample """
    names = [n.rstrip('/') for n in list_html_dir(sample_url(root, sample) + '/') if n.endswith('/')]
    return sorted(n for n in names if not legacy_benchmark_re.match(n) and n not in skip_dirs)


def walk_files(url, depth=3):
    """ (url, filename) for every file at or below a directory, to a bounded depth.  bounded
    because some benchmark directories contain large supplementary trees """
    try:
        names = list_html_dir(url)
    except Exception as e:
        sys.stderr.write('  WARNING: could not list {}: {}\n'.format(url, e))
        return
    for name in names:
        if name.endswith('/'):
            if depth > 1 and name.rstrip('/') not in skip_dirs:
                yield from walk_files(url + name, depth - 1)
        else:
            yield url, name


def detect_ref(url, filename, benchmark):
    """ which reference a benchmark file is against, from its path, then its name, then the
    override table """
    parts = [urllib.parse.unquote(p) for p in url.split('/')]
    for ref in known_refs:
        if ref in parts:
            return ref
    for ref in known_refs:
        if '_{}_'.format(ref) in filename or '_{}.'.format(ref) in filename:
            return ref
    return benchmark_ref_overrides.get(benchmark)


def truthset_label(filename, sample, ref, benchmark):
    """ a short name for this benchmark, for the report's truthset column.

    Drop the tokens that only restate what other columns already say -- the sample, the reference,
    the "1_22" chromosome range, the word "benchmark" -- and keep the rest.  So
    HG002_GRCh38_1_22_v4.2.1_benchmark becomes v4.2.1, and HG002_GRCh38_v5.0q_smvar becomes
    v5.0q_smvar, which is the distinction that actually matters in the report.
    """
    stem = filename[:-len('.vcf.gz')] if filename.endswith('.vcf.gz') else filename
    noise = set([sample.lower(), str(ref).lower(), '1', '22', 'benchmark', 'noinconsistent', ''])
    tokens = [t for t in stem.split('_') if t.lower() not in noise]
    label = '_'.join(tokens).strip('_')
    # nothing distinctive left in the filename, so fall back to the directory it came from
    if len(label) < 3:
        label = benchmark
    return label


def pair_bed(vcf_name, bed_names):
    """ the BED belonging to a VCF: the one sharing the longest common prefix.  GIAB is consistent
    about this even though it is not consistent about anything else """
    stem = vcf_name[:-len('.vcf.gz')] if vcf_name.endswith('.vcf.gz') else vcf_name
    best, best_len = None, 0
    for bed in bed_names:
        shared = len(os.path.commonprefix([stem, bed]))
        if shared > best_len:
            best, best_len = bed, shared
    # a couple of characters in common is coincidence, not a pairing
    if best is not None and best_len >= min(8, len(stem)):
        return best
    return None


def scan_sample(root, sample, benchmarks, ref):
    """ one row per (benchmark, reference-matching VCF) for a sample """
    rows = []
    for benchmark in benchmarks:
        base = '{}/{}/'.format(sample_url(root, sample), benchmark)
        by_dir = {}
        for dir_url, name in walk_files(base):
            by_dir.setdefault(dir_url, []).append(name)

        for dir_url, names in sorted(by_dir.items()):
            vcfs = [n for n in names if n.endswith('.vcf.gz')]
            beds = [n for n in names if n.endswith('.bed') or n.endswith('.bed.gz')]
            for vcf in sorted(vcfs):
                file_ref = detect_ref(dir_url, vcf, benchmark)
                if file_ref != ref:
                    continue
                bed = pair_bed(vcf, beds)
                if not bed:
                    sys.stderr.write('  WARNING: {} {} has no BED alongside {}, skipping -- every '
                                     'comparison tool needs confident regions\n'.format(
                                         sample, benchmark, vcf))
                    continue
                rows.append({'sample': sample,
                             'ref': ref,
                             'truthset': truthset_label(vcf, sample, ref, benchmark),
                             'vcf': dir_url + urllib.parse.quote(vcf),
                             'bed': dir_url + urllib.parse.quote(bed)})
    return rows


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--out', '-o', required=True, help='output manifest CSV')
    ap.add_argument('--samples',
                    help='samples to include (comma-separated, or a file with one per line). '
                    'default is all seven: {}'.format(', '.join(sorted(giab_samples))))
    ap.add_argument('--benchmarks', default=default_benchmarks,
                    help='benchmark set(s) to include, comma separated, or "all" for every '
                    'non-superseded set each sample publishes [{}]'.format(default_benchmarks))
    ap.add_argument('--ref', default=default_ref,
                    help='reference the benchmark is against, one of {} [{}]'.format(
                        ', '.join(known_refs), default_ref))
    ap.add_argument('--root', default=default_root, help='GIAB release tree [{}]'.format(default_root))
    ap.add_argument('--download', metavar='DIR',
                    help='download everything to DIR and write local paths into the manifest')
    ap.add_argument('--threads', type=int, default=8, help='parallel downloads [8]')
    ap.add_argument('--list-benchmarks', action='store_true',
                    help='just print what each sample publishes, and exit')
    args = ap.parse_args()

    samples = parse_sample_list(args.samples) or sorted(giab_samples)
    unknown = [s for s in samples if s not in giab_samples]
    if unknown:
        raise SystemExit('not GIAB benchmark samples: {}. Known samples are: {}'.format(
            ', '.join(sorted(unknown)), ', '.join(sorted(giab_samples))))
    samples = sorted(samples)

    if args.list_benchmarks:
        for sample in samples:
            sys.stderr.write('{}: {}\n'.format(sample, ' '.join(list_benchmarks(args.root, sample))))
        return

    rows = []
    for sample in samples:
        if args.benchmarks == 'all':
            benchmarks = list_benchmarks(args.root, sample)
        else:
            benchmarks = [b.strip() for b in args.benchmarks.split(',') if b.strip()]
        sys.stderr.write('Scanning {} ({})\n'.format(sample, ', '.join(benchmarks)))
        sample_rows = scan_sample(args.root, sample, benchmarks, args.ref)
        if not sample_rows:
            sys.stderr.write('  WARNING: no {} benchmark found for {}\n'.format(args.ref, sample))
        rows += sample_rows

    if not rows:
        raise SystemExit('found nothing for reference {} under {}. Check --ref (the SV Tier1 set '
                         'is GRCh37 only) and --benchmarks.'.format(args.ref, args.root))

    # two benchmarks can produce the same label for one sample only if the tree changes shape
    # under us; catch it here rather than as a confusing duplicate in the report
    seen = set()
    for row in rows:
        key = (row['sample'], row['truthset'])
        if key in seen:
            raise SystemExit('derived the same truthset name "{}" twice for {}. Narrow --benchmarks '
                             'or report this.'.format(row['truthset'], row['sample']))
        seen.add(key)

    if args.download:
        rows = download_rows(rows, args.download, args.threads)

    write_manifest(rows, args.out)
    sys.stderr.write('Wrote {} truth set(s) over {} sample(s) to {}\n'.format(
        len(rows), len(set(r['sample'] for r in rows)), args.out))
    for sample in sorted(set(r['sample'] for r in rows)):
        labels = [r['truthset'] for r in rows if r['sample'] == sample]
        sys.stderr.write('  {}: {}\n'.format(sample, ', '.join(sorted(labels))))


if __name__ == '__main__':
    main()
