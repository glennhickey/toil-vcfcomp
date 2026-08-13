#!/usr/bin/env python3
"""
Pull the headline precision/recall/F1 numbers out of each comparison tool's own summary file and
normalise them into one schema, so that a few hundred per-sample runs can be read as one table.

vcfcomp's *-breakdown functions give per-chromosome FP/FN counts, which is the wrong shape for
comparing across samples.  Each tool does report the headline metrics, but every one of them in a
different format:

  aardvark   summary.tsv     tab separated, already has recall/precision/f1 per variant_type
  vcfeval    summary.txt     fixed-width table, one row per threshold; we want the "None" row
  truvari    summary.json    json (older truvari writes the same json to summary.txt)
  hap.py     *.summary.csv   csv, one row per (Type, Filter)

Aardvark's column set is the superset, so everything else is mapped onto it.
"""

import os, sys
import csv
import glob
import json
import statistics

# aardvark's columns, with the compare_label replaced by sample/tool
report_header = ['callset', 'truthset', 'sample', 'tool', 'comparison', 'region_label', 'filter', 'variant_type',
                 'truth_total', 'truth_tp', 'truth_fn', 'query_total', 'query_tp', 'query_fp',
                 'metric_recall', 'metric_precision', 'metric_f1', 'gt_concordance']

chrom_header = ['callset', 'truthset', 'sample', 'tool', 'metric', 'contig', 'count']

metric_columns = ['metric_recall', 'metric_precision', 'metric_f1']


def _row(callset, truthset, sample, tool, **kwargs):
    row = dict((column, '') for column in report_header)
    row['callset'] = callset
    row['truthset'] = truthset
    row['sample'] = sample
    row['tool'] = tool
    row.update(kwargs)
    return row


def _num(value):
    """ '' and 'nan' both mean "the tool did not report this" """
    try:
        value = float(value)
    except (TypeError, ValueError):
        return None
    return None if value != value else value


def parse_aardvark(out_dir, callset, truthset, sample):
    summary = os.path.join(out_dir, 'summary.tsv')
    if not os.path.isfile(summary):
        return []
    rows = []
    with open(summary, newline='') as summary_file:
        for record in csv.DictReader(summary_file, delimiter='\t'):
            rows.append(_row(callset, truthset, sample, 'aardvark',
                             comparison=record.get('comparison', ''),
                             region_label=record.get('region_label', ''),
                             filter=record.get('filter', ''),
                             variant_type=record.get('variant_type', ''),
                             truth_total=record.get('truth_total', ''),
                             truth_tp=record.get('truth_tp', ''),
                             truth_fn=record.get('truth_fn', ''),
                             query_total=record.get('query_total', ''),
                             query_tp=record.get('query_tp', ''),
                             query_fp=record.get('query_fp', ''),
                             metric_recall=record.get('metric_recall', ''),
                             metric_precision=record.get('metric_precision', ''),
                             metric_f1=record.get('metric_f1', '')))
    return rows


def parse_vcfeval(out_dir, callset, truthset, sample):
    summary = os.path.join(out_dir, 'summary.txt')
    if not os.path.isfile(summary):
        return []
    with open(summary) as summary_file:
        lines = [line for line in summary_file if line.strip() and not line.startswith('---')]
    if len(lines) < 2:
        return []
    # rows are Threshold, True-pos-baseline, True-pos-call, False-pos, False-neg, Precision,
    # Sensitivity, F-measure.  vcfeval also emits a best-threshold row, but the unfiltered "None"
    # row is the one that corresponds to what the other tools report
    data = [line.split() for line in lines[1:]]
    fields = None
    for record in data:
        if record and record[0] == 'None':
            fields = record
    if fields is None:
        fields = data[-1]
    if len(fields) < 8:
        return []

    truth_tp, query_tp, query_fp, truth_fn = (int(f) for f in fields[1:5])
    return [_row(callset, truthset, sample, 'vcfeval',
                 comparison='GT', region_label='ALL', filter='ALL', variant_type='ALL',
                 truth_total=truth_tp + truth_fn,
                 truth_tp=truth_tp,
                 truth_fn=truth_fn,
                 query_total=query_tp + query_fp,
                 query_tp=query_tp,
                 query_fp=query_fp,
                 metric_precision=fields[5],
                 metric_recall=fields[6],
                 metric_f1=fields[7])]


def parse_truvari(out_dir, callset, truthset, sample):
    # truvari v5 writes summary.json; older versions put the same json in summary.txt
    summary = None
    for name in ['summary.json', 'summary.txt']:
        if os.path.isfile(os.path.join(out_dir, name)):
            summary = os.path.join(out_dir, name)
            break
    if not summary:
        return []
    with open(summary) as summary_file:
        record = json.load(summary_file)

    def get(*names):
        # truvari renamed TP-call/call cnt to TP-comp/comp cnt between v3 and v5.  accept both,
        # and complain rather than silently emitting a blank column if it renames them again
        for name in names:
            if name in record:
                value = record[name]
                return '' if value is None else value
        sys.stderr.write('  WARNING: truvari summary has none of {} -- the corresponding report '
                         'column will be blank\n'.format('/'.join(names)))
        return ''

    # NOT 'GT': truvari does not require a genotype match for a TP, it reports genotype agreement
    # separately as gt_concordance.  labelling these rows 'GT' like aardvark's genotype-aware
    # comparison would make a filtered side-by-side read as if truvari simply did better
    return [_row(callset, truthset, sample, 'truvari',
                 comparison='ALLELE', region_label='ALL', filter='ALL', variant_type='SV',
                 truth_total=get('base cnt'),
                 truth_tp=get('TP-base'),
                 truth_fn=get('FN'),
                 query_total=get('comp cnt', 'call cnt'),
                 query_tp=get('TP-comp', 'TP-call'),
                 query_fp=get('FP'),
                 metric_recall=get('recall'),
                 metric_precision=get('precision'),
                 metric_f1=get('f1'),
                 gt_concordance=get('gt_concordance'))]


def parse_happy(out_dir, callset, truthset, sample):
    # hap.py names its outputs after the basename of the output directory
    matches = glob.glob(os.path.join(out_dir, '*.summary.csv'))
    if not matches:
        return []
    rows = []
    with open(matches[0], newline='') as summary_file:
        for record in csv.DictReader(summary_file):
            raw_total = _num(record.get('QUERY.TOTAL'))
            query_fp = _num(record.get('QUERY.FP'))
            query_unk = _num(record.get('QUERY.UNK'))
            # hap.py does not report a query TP directly, and UNK (outside the confident regions)
            # is neither a TP nor an FP, so it has to come out of the total.  QUERY.TOTAL itself
            # includes UNK, so it is subtracted there too -- otherwise query_total would not equal
            # query_tp + query_fp the way it does for the other three tools, and precision
            # recomputed from this table would come out badly wrong for hap.py rows only
            query_tp = ''
            query_total = record.get('QUERY.TOTAL', '')
            if None not in (raw_total, query_fp, query_unk):
                query_tp = int(raw_total - query_fp - query_unk)
                query_total = int(raw_total - query_unk)
            rows.append(_row(callset, truthset, sample, 'happy',
                             comparison='GT', region_label='ALL',
                             filter=record.get('Filter', ''),
                             variant_type=record.get('Type', ''),
                             truth_total=record.get('TRUTH.TOTAL', ''),
                             truth_tp=record.get('TRUTH.TP', ''),
                             truth_fn=record.get('TRUTH.FN', ''),
                             query_total=query_total,
                             query_tp=query_tp,
                             query_fp=record.get('QUERY.FP', ''),
                             metric_recall=record.get('METRIC.Recall', ''),
                             metric_precision=record.get('METRIC.Precision', ''),
                             metric_f1=record.get('METRIC.F1_Score', '')))
    return rows


summary_parsers = {'aardvark': parse_aardvark,
                   'vcfeval': parse_vcfeval,
                   'truvari': parse_truvari,
                   'happy': parse_happy}


def parse_summary(tool, out_dir, callset, truthset, sample):
    """ headline metrics for one tool's output directory, as a list of report_header dicts """
    if tool not in summary_parsers:
        raise RuntimeError('no summary parser for tool "{}"'.format(tool))
    return summary_parsers[tool](out_dir, callset, truthset, sample)


def chromsplit_to_long(table, callset, truthset, sample, tool):
    """ vcfcomp's *_chromsplit tables are wide (one column per contig) and the contig set differs
    between samples, so they cannot just be stacked.  melt them into long form instead """
    if not table or len(table) < 2:
        return []
    contigs = table[0][1:]
    rows = []
    for record in table[1:]:
        metric = record[0]
        for contig, count in zip(contigs, record[1:]):
            rows.append({'callset': callset, 'truthset': truthset, 'sample': sample,
                         'tool': tool, 'metric': metric,
                         'contig': contig, 'count': count})
    return rows


def write_tsv(rows, header, path):
    with open(path, 'w', newline='') as out_file:
        writer = csv.DictWriter(out_file, fieldnames=header, delimiter='\t', extrasaction='ignore')
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def summarize(rows):
    """ collapse the per-sample rows into one row per (tool, comparison, region, filter, type) """
    groups = {}
    for row in rows:
        # callset is part of the key: the whole point of running several is to compare them, so
        # averaging across them would destroy the answer
        key = (row['callset'], row['truthset'], row['tool'], row['comparison'], row['region_label'],
               row['filter'], row['variant_type'])
        groups.setdefault(key, []).append(row)

    header = ['callset', 'truthset', 'tool', 'comparison', 'region_label', 'filter', 'variant_type',
              'samples']
    for metric in metric_columns:
        name = metric.replace('metric_', '')
        header += ['{}_{}'.format(name, stat) for stat in ['mean', 'median', 'min', 'max']]

    out_rows = []
    for key in sorted(groups):
        group = groups[key]
        row = dict(zip(header, list(key) + [len(group)]))
        for metric in metric_columns:
            name = metric.replace('metric_', '')
            values = [v for v in (_num(r[metric]) for r in group) if v is not None]
            if values:
                row['{}_mean'.format(name)] = '{:.6f}'.format(statistics.fmean(values))
                row['{}_median'.format(name)] = '{:.6f}'.format(statistics.median(values))
                row['{}_min'.format(name)] = '{:.6f}'.format(min(values))
                row['{}_max'.format(name)] = '{:.6f}'.format(max(values))
            else:
                for stat in ['mean', 'median', 'min', 'max']:
                    row['{}_{}'.format(name, stat)] = ''
        out_rows.append(row)
    return header, out_rows


# each tool reports one comparison that means "did we get the call right".  they do not agree on
# what to call it, and mixing them in one table would compare unlike things
primary_comparison = {'aardvark': 'GT', 'vcfeval': 'GT', 'happy': 'GT', 'truvari': 'ALLELE'}


def headline_rows(summary_rows, comparison=None):
    """ the subset of the summary worth putting in front of someone: the primary comparison for
    each tool, unstratified by region or filter """
    keep = []
    for row in summary_rows:
        if row['region_label'] not in ('ALL', ''):
            continue
        if row['filter'] not in ('ALL', ''):
            continue
        wanted = comparison or primary_comparison.get(row['tool'], 'GT')
        if row['comparison'] != wanted:
            continue
        keep.append(row)
    return keep


def _variant_type_order(names):
    """ ALL first, then alphabetical -- the overall number is the one people want first """
    return ['ALL'] + sorted(n for n in names if n != 'ALL')


def build_comparison(summary_rows, callset_order, metric='f1', comparison=None):
    """ head-to-head: one row per (tool, truthset, variant type), one column per callset, plus the
    delta against the first callset and which callset won """
    rows = headline_rows(summary_rows, comparison)
    column = '{}_mean'.format(metric)

    cells = {}
    callsets = []
    for row in rows:
        if row['callset'] not in callsets:
            callsets.append(row['callset'])
        cells[(row['tool'], row['truthset'], row['variant_type'], row['callset'])] = _num(row[column])

    # the order the callsets were given on the command line, so the baseline is the first --queryVcf
    if callset_order:
        callsets = [c for c in callset_order if c in callsets] + \
                   [c for c in callsets if c not in callset_order]
    else:
        callsets = sorted(callsets)

    keys = sorted(set((t, ts) for (t, ts, _v, _c) in cells))
    out = []
    for tool, truthset in keys:
        variant_types = set(v for (t, ts, v, _c) in cells if t == tool and ts == truthset)
        for variant_type in _variant_type_order(variant_types):
            values = [cells.get((tool, truthset, variant_type, c)) for c in callsets]
            present = [(c, v) for c, v in zip(callsets, values) if v is not None]
            if not present:
                continue
            row = {'tool': tool, 'truthset': truthset, 'variant_type': variant_type,
                   'metric': metric}
            for callset, value in zip(callsets, values):
                row[callset] = '' if value is None else '{:.6f}'.format(value)
            baseline = values[0]
            row['delta'] = ''
            if baseline is not None and len(present) > 1 and values[-1] is not None:
                row['delta'] = '{:+.6f}'.format(values[-1] - baseline)
            row['best'] = max(present, key=lambda cv: cv[1])[0] if len(present) > 1 else ''
            out.append(row)

    header = ['tool', 'truthset', 'variant_type', 'metric'] + callsets + ['delta', 'best']
    return header, out, callsets


def _fmt_table(header, rows, aligns=None):
    """ a whitespace-aligned block, since this file is meant to be read rather than parsed """
    widths = [len(h) for h in header]
    for row in rows:
        for i, cell in enumerate(row):
            widths[i] = max(widths[i], len(str(cell)))
    aligns = aligns or ['<'] * len(header)
    out = ['  '.join('{:{}{}}'.format(h, aligns[i], widths[i]) for i, h in enumerate(header))]
    for row in rows:
        out.append('  '.join('{:{}{}}'.format(str(c), aligns[i], widths[i])
                             for i, c in enumerate(row)))
    return out


def write_summary_text(rows, summary_rows, path, callset_order=None, comparison=None):
    """ the file to open first: what ran, the headline numbers, and which callset won """
    lines = []
    heads = headline_rows(summary_rows, comparison)
    callsets = callset_order or sorted(set(r['callset'] for r in rows))
    truthsets = sorted(set(r['truthset'] for r in rows))
    tools = sorted(set(r['tool'] for r in rows))
    n_samples = len(set(r['sample'] for r in rows))

    lines.append('toil-vcfcomp summary')
    lines.append('====================')
    lines.append('')
    lines.append('callsets  : {}{}'.format(', '.join(callsets),
                                           '   (baseline: {})'.format(callsets[0]) if len(callsets) > 1 else ''))
    lines.append('truthsets : {}'.format(', '.join(truthsets)))
    lines.append('samples   : {}'.format(n_samples))
    lines.append('tools     : {}'.format(', '.join('{} (comparison={})'.format(
        t, comparison or primary_comparison.get(t, 'GT')) for t in tools)))
    lines.append('')

    header, compare_rows, compare_callsets = build_comparison(summary_rows, callsets,
                                                              comparison=comparison)

    if len(compare_callsets) > 1:
        wins = {}
        for row in compare_rows:
            if row['best']:
                wins[row['best']] = wins.get(row['best'], 0) + 1
        total = sum(wins.values())
        if total:
            winner = max(wins, key=lambda c: wins[c])
            lines.append('VERDICT: {} has the better F1 on {} of {} (tool, truthset, variant type) '
                         'combinations.'.format(winner, wins[winner], total))
            others = ', '.join('{}: {}'.format(c, wins.get(c, 0)) for c in compare_callsets if c != winner)
            if others:
                lines.append('         ({})'.format(others))
        lines.append('')

    # headline: the ALL row, with recall and precision alongside F1 so a precision/recall swap
    # is visible rather than hidden inside the F1
    table = []
    for row in sorted(heads, key=lambda r: (r['tool'], r['truthset'],
                                            callsets.index(r['callset']) if r['callset'] in callsets else 0)):
        if row['variant_type'] != 'ALL':
            continue
        table.append([row['tool'], row['truthset'], row['callset'], row['samples'],
                      row['recall_mean'][:7], row['precision_mean'][:7], row['f1_mean'][:7]])
    if table:
        lines.append('Headline  (variant_type=ALL, mean across samples)')
        lines += ['  ' + l for l in _fmt_table(
            ['tool', 'truthset', 'callset', 'n', 'recall', 'precis.', 'F1'], table)]
        lines.append('')

    if compare_rows:
        title = 'F1 by variant type' + (' (delta = last callset - baseline)' if len(compare_callsets) > 1 else '')
        lines.append(title)
        display = ['truthset', 'variant_type'] + compare_callsets + (
            ['delta', 'best'] if len(compare_callsets) > 1 else [])
        table = []
        for row in compare_rows:
            cells = [row['truthset'], row['variant_type']] + \
                    [row[c][:7] if row[c] else '-' for c in compare_callsets]
            if len(compare_callsets) > 1:
                cells += [row['delta'][:9] if row['delta'] else '-', row['best'] or '-']
            table.append(cells)
        lines += ['  ' + l for l in _fmt_table(display, table)]
        lines.append('')

    lines.append('Full tables: report.tsv (per sample), report.summary.tsv (across samples), '
                 'report.compare.tsv (head to head), report.chrom.tsv (per contig).')

    with open(path, 'w') as out_file:
        out_file.write('\n'.join(lines) + '\n')
    return header, compare_rows


def write_report(rows, chrom_rows, out_dir, prefix='report', callset_order=None, comparison=None):
    """ write the report tables and return their paths """
    paths = {}

    rows = sorted(rows, key=lambda r: (r['callset'], r['truthset'], r['sample'], r['tool'],
                                       r['variant_type'], r['filter']))
    paths['report'] = os.path.join(out_dir, prefix + '.tsv')
    write_tsv(rows, report_header, paths['report'])

    header, summary_rows = summarize(rows)
    paths['summary'] = os.path.join(out_dir, prefix + '.summary.tsv')
    write_tsv(summary_rows, header, paths['summary'])

    paths['text'] = os.path.join(out_dir, prefix + '.txt')
    compare_header, compare_rows = write_summary_text(rows, summary_rows, paths['text'],
                                                      callset_order=callset_order,
                                                      comparison=comparison)
    if compare_rows:
        paths['compare'] = os.path.join(out_dir, prefix + '.compare.tsv')
        write_tsv(compare_rows, compare_header, paths['compare'])

    if chrom_rows:
        paths['chrom'] = os.path.join(out_dir, prefix + '.chrom.tsv')
        write_tsv(sorted(chrom_rows, key=lambda r: (r['callset'], r['truthset'], r['sample'],
                                                    r['tool'], r['metric'], r['contig'])),
                  chrom_header, paths['chrom'])

    return paths


def main():
    """ standalone: point it at a directory of per-sample tool outputs and rebuild the report """
    import argparse
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--dir', required=True,
                    help='directory of <callset>/<truthset>/<sample>/<tool> output directories, as '
                    'written by toil-vcfcomp. shallower layouts from earlier versions also work')
    ap.add_argument('--out-dir', help='where to write the report tables [--dir]')
    ap.add_argument('--prefix', default='report', help='output file prefix [report]')
    args = ap.parse_args()

    # find any directory named after a tool, then read sample/truthset/callset off its parents.
    # walking up like this means every layout this tool has ever written still parses
    root = os.path.abspath(args.dir)

    def parent_name(path, levels):
        for _ in range(levels):
            if os.path.abspath(path) == root:
                return '', path
            path = os.path.dirname(path)
        return ('' if os.path.abspath(path) == root else os.path.basename(path)), path

    rows = []
    for dirpath, dirnames, _ in os.walk(root):
        for tool in sorted(d for d in dirnames if d in summary_parsers):
            sample = os.path.basename(dirpath)
            truthset, up = parent_name(dirpath, 1)
            callset, _up2 = parent_name(dirpath, 2)
            rows += parse_summary(tool, os.path.join(dirpath, tool), callset, truthset, sample)

    if not rows:
        raise SystemExit('no tool summaries found under {}'.format(args.dir))

    paths = write_report(rows, [], args.out_dir if args.out_dir else args.dir, args.prefix)
    for name in sorted(paths):
        sys.stderr.write('Wrote {}\n'.format(paths[name]))


if __name__ == '__main__':
    main()
