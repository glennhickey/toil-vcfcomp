#!/usr/bin/env python3
"""
The truth manifest shared by the scrapers and the workflow, plus the downloading that goes with it.

A manifest is a CSV with one row per truth set:

  sample,ref,truthset,vcf,bed

`truthset` names which benchmark the row is (dipcall, v4.2.1, CMRG_smallvar_v1.00, ...).  It is a
reporting dimension: one run can compare several call sets against several truth sets at once.
Manifests written before that column existed used `caller` instead, and are still accepted.
"""

import os, sys
import csv
import re
import shutil
import subprocess
import time
import urllib.request
import urllib.parse
from concurrent.futures import ThreadPoolExecutor

manifest_header = ['sample', 'ref', 'truthset', 'vcf', 'bed']
# what the truthset column used to be called
legacy_label_columns = ['caller']


def list_html_dir(url, timeout=90):
    """ names listed by an apache/nginx style directory index.  the GIAB FTP tree is served this
    way and has no machine-readable index, so this is the only way to enumerate it """
    with urllib.request.urlopen(url, timeout=timeout) as response:
        body = response.read().decode('utf-8', 'replace')
    names = []
    for href in re.findall(r'href="([^"]+)"', body):
        if href.startswith(('/', '?', 'http')) or href.startswith('..'):
            continue
        names.append(urllib.parse.unquote(href))
    return names


def write_manifest(rows, path):
    with open(path, 'w', newline='') as manifest_file:
        writer = csv.DictWriter(manifest_file, fieldnames=manifest_header)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)

def read_manifest(path):
    """ read a manifest back.  shared with the toil workflow """
    with open(path, newline='') as manifest_file:
        reader = csv.DictReader(manifest_file)
        fieldnames = reader.fieldnames or []
        # accept a manifest written before the truthset column existed
        label_column = 'truthset'
        if 'truthset' not in fieldnames:
            for legacy in legacy_label_columns:
                if legacy in fieldnames:
                    label_column = legacy
                    break
        required = [c if c != 'truthset' else label_column for c in manifest_header]
        missing = [c for c in required if c not in fieldnames]
        if missing:
            raise SystemExit('{} is missing required column(s): {}. Expected a header of: {}'.format(
                path, ', '.join(missing), ','.join(manifest_header)))
        rows = []
        for row in reader:
            row = dict(row)
            row['truthset'] = row.get('truthset') or row.get(label_column) or 'truth'
            rows.append(row)

    if not rows:
        raise SystemExit('{} contains no rows'.format(path))

    # (sample, truthset) is the key the workflow fans out on, so a duplicate would silently run
    # the same comparison twice and double-count that sample in the summary
    seen = {}
    for row in rows:
        key = (row['sample'], row['truthset'])
        if key in seen:
            raise SystemExit('{} has duplicate rows for sample {} / truthset {}. Each pair must '
                             'appear once.'.format(path, key[0], key[1]))
        seen[key] = True
    return rows

def remote_size(url):
    """ Content-Length for a url, or None if the server does not say """
    if url.startswith('s3://'):
        return None
    try:
        request = urllib.request.Request(url, method='HEAD')
        with urllib.request.urlopen(request, timeout=60) as response:
            length = response.headers.get('Content-Length')
            return int(length) if length is not None else None
    except Exception:
        return None

def download_url(url, path, retries=5):
    """ fetch one url to a local path, atomically.

    Two things are being defended against.  A half-written file that still parses is worse than no
    file at all, so the transfer lands on .part and is renamed only once complete.  And a body that
    ends early after a clean close does NOT raise -- the sized read just returns b'' and
    copyfileobj stops -- so the transfer is only trusted when the bytes written match
    Content-Length.  Raises RuntimeError rather than SystemExit: Toil special-cases SystemExit and
    skips its traceback logging, so the failing URL would never reach the job log.
    """
    if url.startswith('s3://'):
        # the AWS CLI, which is not a dependency of this package.  toil-vcfcomp preflights for it
        if not shutil.which('aws'):
            raise RuntimeError('{} is an s3:// url, which needs the AWS CLI ("aws") on the PATH. '
                               'It is not installed. Rebuild the manifest without --s3-urls to get '
                               'https urls, which need no extra tooling.'.format(url))
        subprocess.run(['aws', 's3', 'cp', '--no-sign-request', '--quiet', url, path], check=True)
        return path

    for attempt in range(retries):
        try:
            written = 0
            with urllib.request.urlopen(url, timeout=300) as response, open(path + '.part', 'wb') as out_file:
                # the GET's own Content-Length, so no extra HEAD round trip per file and the check
                # still works against a server that does not answer HEAD
                declared = response.headers.get('Content-Length')
                expected = int(declared) if declared is not None else None
                while True:
                    chunk = response.read(1 << 20)
                    if not chunk:
                        break
                    out_file.write(chunk)
                    written += len(chunk)
            if expected is not None and written != expected:
                raise IOError('short read: got {} bytes, expected {}'.format(written, expected))
            os.rename(path + '.part', path)
            return path
        except Exception as e:
            try:
                os.remove(path + '.part')
            except OSError:
                pass
            if attempt == retries - 1:
                raise RuntimeError('failed to download {} after {} attempts: {}'.format(url, retries, e))
            sys.stderr.write('  retrying {} after error: {}\n'.format(url, e))
            time.sleep(min(2 ** attempt, 30))

def is_complete(url, local):
    """ a previously downloaded file is only reusable if it is the whole file.  "it exists" is not
    good enough: a truncated leftover from an earlier run would be silently trusted forever """
    if not os.path.isfile(local):
        return False
    expected = remote_size(url)
    if expected is None:
        return True
    if os.path.getsize(local) == expected:
        return True
    sys.stderr.write('  {} is {} bytes but the server says {} -- redownloading\n'.format(
        os.path.basename(local), os.path.getsize(local), expected))
    return False

def download_rows(rows, out_dir, threads):
    """ download every vcf (+ its index) and bed, and rewrite the rows to point at the copies """
    os.makedirs(out_dir, exist_ok=True)
    jobs = []
    for row in rows:
        for column in ['vcf', 'bed']:
            if not row[column]:
                continue
            urls = [row[column]] + ([row[column] + '.tbi'] if column == 'vcf' else [])
            for url in urls:
                local = os.path.join(out_dir, os.path.basename(urllib.parse.urlparse(url).path))
                jobs.append((url, local))

    sys.stderr.write('Downloading {} files to {}\n'.format(len(jobs), out_dir))
    with ThreadPoolExecutor(max_workers=threads) as pool:
        futures = []
        for url, local in jobs:
            if is_complete(url, local):
                continue
            futures.append(pool.submit(download_url, url, local))
        for future in futures:
            future.result()

    for row in rows:
        for column in ['vcf', 'bed']:
            if row[column]:
                local = os.path.join(out_dir, os.path.basename(urllib.parse.urlparse(row[column]).path))
                row[column] = os.path.abspath(local)
    return rows

def parse_sample_list(value):
    """ accept either a comma-separated list or a file with one sample per line """
    if not value:
        return None
    if os.path.isfile(value):
        with open(value) as sample_file:
            return set(line.strip() for line in sample_file if line.strip())
    return set(s.strip() for s in value.split(',') if s.strip())

