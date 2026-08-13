#!/usr/bin/env python3
"""
Turn an HPRC variant-call submission on s3 into a truth set for toil-vcfcomp.

The Yale R2 submission stores one directory per sample, laid out as

  samples/<SAMPLE>/variant_calls/<REF>/<CALLER>/<SAMPLE>.<REF>.<CALLER>.vcf.gz

with dipcall additionally providing <SAMPLE>.<REF>.dipcall.bed, the confident regions the
comparison has to be restricted to.  This scrapes that layout and writes a manifest CSV in the
same shape as the HPRC annotation index files:

  sample,ref,caller,vcf,bed

The manifest is what toil-vcfcomp consumes.  --merge additionally builds one big multi-sample
truth VCF, which is a useful artifact on its own but is deliberately NOT the workflow's input:
a merged VCF cannot carry the per-sample confident BEDs, and re-subsetting one sample out of it
is slower than just reading that sample's original file.

URLs are written as https by default, because they are fetched with the standard library and need
neither credentials nor extra tooling.  The s3:// form (--s3-urls) is fetched by shelling out to
the AWS CLI, which is not a dependency of this package and is routinely absent from compute nodes.

Usage:

  hprc-truthset --ref GRCh38_no_alt --caller dipcall -o truth.index.csv

  hprc-truthset --ref GRCh38_no_alt --download ./truth -o truth.index.csv

  hprc-truthset --ref GRCh38_no_alt --download ./truth -o truth.index.csv \
                --merge truth.merged.vcf.gz --threads 16

Needs bcftools and tabix on PATH for --merge.  Python 3, stdlib only.
"""

import os, sys
import argparse
import csv
import subprocess
import shutil
import time
import urllib.request
import urllib.parse
import xml.etree.ElementTree as ET
from concurrent.futures import ThreadPoolExecutor

default_bucket = 'human-pangenomics'
default_region = 'us-west-2'
default_prefix = 'submissions/40399FDD-59DE-43D1-B3A3-DFF0C6E64FAC--YALE_VARIANT_CALLS_R2/'
default_ref = 'GRCh38_no_alt'
default_caller = 'dipcall'

s3_ns = '{http://s3.amazonaws.com/doc/2006-03-01/}'

from toil_vcfcomp.manifest import (manifest_header, read_manifest, write_manifest,
                                   download_rows, parse_sample_list)

# --- s3 ---

def bucket_endpoint(bucket, region):
    return 'https://{}.s3.{}.amazonaws.com'.format(bucket, region)

def key_url(bucket, region, key, s3_urls=False):
    if s3_urls:
        return 's3://{}/{}'.format(bucket, key)
    return '{}/{}'.format(bucket_endpoint(bucket, region), urllib.parse.quote(key))

def list_keys(bucket, region, prefix):
    """ list every key under a prefix using the anonymous s3 REST api.  no credentials and no
    aws cli needed, which matters because neither is a given on a compute node """
    token = None
    n_pages = 0
    while True:
        query = {'list-type': '2', 'prefix': prefix, 'max-keys': '1000'}
        if token:
            query['continuation-token'] = token
        url = '{}/?{}'.format(bucket_endpoint(bucket, region), urllib.parse.urlencode(query))
        with urllib.request.urlopen(url, timeout=120) as response:
            root = ET.fromstring(response.read())
        for contents in root.findall(s3_ns + 'Contents'):
            yield contents.find(s3_ns + 'Key').text
        n_pages += 1
        sys.stderr.write('  listed {} pages of {}\r'.format(n_pages, prefix))
        truncated = root.find(s3_ns + 'IsTruncated')
        if truncated is None or truncated.text != 'true':
            sys.stderr.write('\n')
            return
        token = root.find(s3_ns + 'NextContinuationToken').text

# --- manifest ---

def scan_submission(bucket, region, prefix, ref, caller, bed_caller, s3_urls):
    """ walk the submission and pick out one vcf (and its confident-region bed) per sample """
    if not prefix.endswith('/'):
        prefix += '/'
    samples_prefix = prefix + 'samples/'

    vcfs = {}
    beds = {}
    for key in list_keys(bucket, region, samples_prefix):
        rest = key[len(samples_prefix):]
        parts = rest.split('/')
        # <SAMPLE>/variant_calls/<REF>/<CALLER>/<FILE>
        if len(parts) != 5 or parts[1] != 'variant_calls' or parts[2] != ref:
            continue
        sample, key_caller, name = parts[0], parts[3], parts[4]
        if key_caller == caller and name.endswith('.vcf.gz'):
            # some callers publish several vcfs per directory (eg svdss ships a .genotyped one
            # alongside the plain one).  take the one named exactly for the caller
            if name == '{}.{}.{}.vcf.gz'.format(sample, ref, caller):
                vcfs[sample] = key_url(bucket, region, key, s3_urls)
        if key_caller == bed_caller and name.endswith('.bed'):
            beds[sample] = key_url(bucket, region, key, s3_urls)

    rows = []
    for sample in sorted(vcfs):
        rows.append({'sample': sample, 'ref': ref, 'truthset': caller,
                     'vcf': vcfs[sample], 'bed': beds.get(sample, '')})
    return rows

# --- merge ---

def vcf_contigs(vcf_path):
    """ contigs that actually have records, in index order """
    result = subprocess.run(['bcftools', 'index', '-s', vcf_path], capture_output=True, text=True, check=True)
    return [line.split('\t')[0] for line in result.stdout.strip().split('\n') if line.strip()]

def merge_truth(rows, out_vcf, work_dir, threads, merge_options):
    """ build one multi-sample truth vcf by merging per chromosome and concatenating.  a single
    whole-genome bcftools merge over a couple of hundred files is needlessly serial """
    vcfs = [row['vcf'] for row in rows]
    for vcf in vcfs:
        if not os.path.isfile(vcf):
            raise SystemExit('--merge needs local files: {} is not on disk. Run with --download first.'.format(vcf))

    os.makedirs(work_dir, exist_ok=True)

    # union of contigs across samples, keeping the first sample's order and appending anything
    # only the others have (chrY shows up in males only)
    contigs = []
    seen = set()
    for vcf in vcfs:
        for contig in vcf_contigs(vcf):
            if contig not in seen:
                seen.add(contig)
                contigs.append(contig)
    sys.stderr.write('Merging {} samples over {} contigs\n'.format(len(vcfs), len(contigs)))

    file_list = os.path.join(work_dir, 'merge-inputs.txt')
    with open(file_list, 'w') as list_file:
        for vcf in vcfs:
            list_file.write(vcf + '\n')

    def merge_contig(contig):
        chunk = os.path.join(work_dir, 'merged.{}.vcf.gz'.format(contig))
        cmd = ['bcftools', 'merge', '-r', contig, '-l', file_list, '-Oz', '-o', chunk]
        if merge_options:
            cmd += merge_options.split()
        sys.stderr.write('  merging {}\n'.format(contig))
        subprocess.run(cmd, check=True)
        subprocess.run(['tabix', '-fp', 'vcf', chunk], check=True)
        return chunk

    with ThreadPoolExecutor(max_workers=threads) as pool:
        chunks = list(pool.map(merge_contig, contigs))

    sys.stderr.write('Concatenating {} chunks into {}\n'.format(len(chunks), out_vcf))
    subprocess.run(['bcftools', 'concat', '-Oz', '-o', out_vcf] + chunks, check=True)
    subprocess.run(['tabix', '-fp', 'vcf', out_vcf], check=True)
    for chunk in chunks:
        os.remove(chunk)
        os.remove(chunk + '.tbi')
    os.remove(file_list)

# --- main ---

def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--out', '-o', required=True, help='output manifest CSV')
    ap.add_argument('--ref', default=default_ref,
                    help='reference the calls were made against [{}]'.format(default_ref))
    ap.add_argument('--caller', default=default_caller,
                    help='variant caller subdirectory to take the truth VCFs from [{}]'.format(default_caller))
    ap.add_argument('--bed-caller', default=default_caller,
                    help='caller subdirectory to take the confident-region BED from. only dipcall '
                    'publishes one, so leave this alone even when --caller is something else [{}]'.format(default_caller))
    ap.add_argument('--bucket', default=default_bucket, help='s3 bucket [{}]'.format(default_bucket))
    ap.add_argument('--region', default=default_region, help='s3 region [{}]'.format(default_region))
    ap.add_argument('--prefix', default=default_prefix,
                    help='submission prefix within the bucket [{}]'.format(default_prefix))
    ap.add_argument('--samples', help='restrict to these samples (comma-separated, or a file with one per line)')
    ap.add_argument('--exclude', help='drop these samples (comma-separated, or a file with one per line)')
    ap.add_argument('--s3-urls', action='store_true',
                    help='write s3:// urls instead of https. these are fetched with the AWS CLI, '
                    'which is not a dependency of this package and is often missing from compute '
                    'nodes; https needs no extra tooling and no credentials')
    ap.add_argument('--download', metavar='DIR',
                    help='download everything to DIR and write local paths into the manifest. use this '
                    'when the compute nodes have no outbound network access')
    ap.add_argument('--merge', metavar='OUT.vcf.gz',
                    help='also build a single multi-sample truth VCF. requires --download, since '
                    'bcftools merge cannot read the inputs over https')
    ap.add_argument('--merge-options', default='',
                    help='extra options for bcftools merge (surround in quotes)')
    ap.add_argument('--threads', type=int, default=8, help='parallel downloads / merges [8]')
    args = ap.parse_args()

    # checked before the scan, so a missing flag does not cost 20 seconds of listing first
    if args.merge and not args.download:
        raise SystemExit('--merge needs --download DIR: bcftools merge cannot read the VCFs over https')
    if args.s3_urls and not shutil.which('aws'):
        sys.stderr.write('  WARNING: --s3-urls writes urls that are fetched with the AWS CLI, and '
                         '"aws" is not on your PATH\n')

    sys.stderr.write('Scanning s3://{}/{}\n'.format(args.bucket, args.prefix))
    rows = scan_submission(args.bucket, args.region, args.prefix, args.ref, args.caller,
                           args.bed_caller, args.s3_urls)
    if not rows:
        raise SystemExit('found no {} VCFs for reference {} under s3://{}/{}. Check --ref and --caller.'.format(
            args.caller, args.ref, args.bucket, args.prefix))

    keep = parse_sample_list(args.samples)
    if keep:
        found = set(row['sample'] for row in rows)
        for sample in sorted(keep - found):
            sys.stderr.write('  WARNING: requested sample {} not found in the submission\n'.format(sample))
        rows = [row for row in rows if row['sample'] in keep]
    drop = parse_sample_list(args.exclude)
    if drop:
        rows = [row for row in rows if row['sample'] not in drop]

    no_bed = [row['sample'] for row in rows if not row['bed']]
    if no_bed:
        sys.stderr.write('  WARNING: {} sample(s) have no {} BED, so they have no confident regions '
                         'of their own: {}\n'.format(len(no_bed), args.bed_caller,
                                                     ', '.join(no_bed[:5]) + (' ...' if len(no_bed) > 5 else '')))

    if args.download:
        rows = download_rows(rows, args.download, args.threads)

    write_manifest(rows, args.out)
    sys.stderr.write('Wrote {} samples to {}\n'.format(len(rows), args.out))

    if args.merge:
        merge_truth(rows, args.merge, os.path.join(args.download, 'merge-work'),
                    args.threads, args.merge_options)
        sys.stderr.write('Wrote merged truth VCF to {}\n'.format(args.merge))

if __name__ == '__main__':
    main()
