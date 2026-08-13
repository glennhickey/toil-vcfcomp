# toil-vcfcomp

Benchmark every sample of a multi-sample VCF against its own truth set, in parallel with
[Toil](https://github.com/DataBiosphere/toil).

`vcfcomp` compares one call set against one truth set. `toil-vcfcomp` runs that comparison once
per sample -- typically a pangenome graph VCF with a few hundred samples in it, against the
per-sample dipcall calls published for HPRC -- and collects the results into a single report.

Two steps:

1. A **truth-set scraper** writes a manifest CSV naming one truth VCF and one confident-region BED
   per row. `hprc-truthset` scrapes an HPRC variant-call submission on s3; `giab-truthset` scrapes
   the GIAB benchmark FTP tree.
2. **`toil-vcfcomp`** takes that manifest plus one or more multi-sample query VCFs, works out which
   samples they have in common, and fans out one comparison job per callset × truthset × sample.

The report has four dimensions, so a single run can answer "which of my graphs is better" and
"how does it look against each benchmark" at the same time:

| dimension | comes from |
|---|---|
| `callset` | each `--queryVcf` |
| `truthset` | each `truthset` value in the manifest |
| `sample` | samples shared between a callset and the truth |
| `tool` | each `--evalType` |

## Installing

Toil is pinned separately so that a Toil release cannot silently change the workflow, and it has
to go in first:

```
python3 -m venv venv && . venv/bin/activate
pip install -U -r toil-requirement.txt
pip install -U .
```

Install into the *same* virtualenv as Toil. Toil launches its workers by running `_toil_worker`
off the `PATH`, so a stale copy earlier on the `PATH` will be picked up instead and every job
will fail with `ModuleNotFoundError: No module named 'toil'`.

Also needed on every worker: `bcftools` (>= 1.10, for the `+split` plugin), `tabix`, `bgzip`,
`samtools` if the reference has no `.fai` yet, and `bedtools` if you pass `--regions`.
The comparison tools themselves come from containers -- see `--binariesMode`.

## Building a truth set

```
hprc-truthset --ref GRCh38_no_alt --caller dipcall -o truth.index.csv
```

That writes one row per sample:

```
sample,ref,truthset,vcf,bed
HG00097,GRCh38_no_alt,dipcall,https://human-pangenomics.s3.../HG00097.GRCh38_no_alt.dipcall.vcf.gz,https://.../HG00097.GRCh38_no_alt.dipcall.bed
```

232 samples, of which **231** overlap the HPRC r2 graph — including **HG002 and HG005**, the two
GIAB benchmark samples that are also in the pangenome. Manifests written before the `truthset`
column existed used `caller` instead; those are still read.

URLs are https rather than `s3://` because they are fetched with the standard library and need
neither credentials nor extra tooling. `--s3-urls` switches, but those are downloaded by shelling
out to the AWS CLI, which is not a dependency of this package and is often absent from compute
nodes. Useful options:

- `--samples` / `--exclude` to subset.
- `--download DIR` to fetch everything up front and write local paths into the manifest. Use this
  when the compute nodes have no outbound network access. Downloads are verified against
  `Content-Length` and staged through a `.part` file, so an interrupted run leaves no
  half-written file that a later run would trust.
- `--merge OUT.vcf.gz` to also build one multi-sample truth VCF. Requires `--download`, since
  `bcftools merge` cannot read its inputs over https. This is a standalone artifact, not the
  workflow's input: a merged VCF cannot carry the per-sample confident BEDs, and pulling one
  sample back out of it is slower than reading that sample's original file.

Both `--ref CHM13v2` and other `--caller` values work, but only dipcall publishes a
confident-region BED. For any other caller, either keep `--bed-caller dipcall` or supply a global
`--regions` BED to the workflow.

## Building a GIAB truth set

```
giab-truthset -o giab.index.csv                             # v4.2.1, all 7 samples, GRCh38
giab-truthset --samples HG002 --benchmarks all -o giab.csv  # every benchmark HG002 has
giab-truthset --list-benchmarks                             # what each sample publishes
```

GIAB publishes far more for HG002 than for anyone else. Each benchmark becomes its own `truthset`
row, so one run scores a call set against all of them side by side:

| sample | truth sets available on GRCh38 |
|---|---|
| HG002 | `v4.2.1`, `CMRG_smallvar_v1.00`, `CMRG_SV_v1.00`, `TandemRepeats_v1.0`, `TandemRepeats_v1.0.1`, `chrXY_smallvar_v1.0`, `v5.0q_smvar`, `v5.0q_stvar`, `MosaicSNVv1.0/v1.1` (+`_GermlineV4.2.1`) |
| HG001, HG003–HG007 | `v4.2.1` |

The GIAB tree is not uniformly shaped — some benchmarks sit directly in the benchmark directory,
some under a per-reference subdirectory, some under a further SmallVariant/StructuralVariant split
— so the scraper walks each directory and pairs a VCF with the BED sharing the longest name prefix,
which is the one thing GIAB is consistent about. A VCF with no BED beside it is skipped with a
warning rather than silently compared against the wrong regions.

Two things to watch:

- **`--ref` matters.** The SV Tier1 v0.6 set is GRCh37-only and simply will not appear under
  `--ref GRCh38`; Q100 (`v5.0q`) is published for GRCh37, GRCh38 and CHM13v2.0.
- **Match the tool to the benchmark.** `v4.2.1`, `CMRG_smallvar` and `chrXY_smallvar` are
  small-variant benchmarks (vcfeval, hap.py, aardvark); `CMRG_SV`, `SVs_Tier1` and `v5.0q_stvar`
  are structural (truvari, aardvark). Running a small-variant tool against an SV benchmark will
  complete and report nonsense — nothing detects that for you.

Only **HG002 and HG005** are in the HPRC r2 graph, so a GIAB run against that graph covers those
two samples; the rest of the manifest is reported as truth-only and skipped.

## Running the comparison

```
toil-vcfcomp ./jobstore \
    --queryVcf hprc-r2-mc-grch38.vcf.gz \
    --truthIndex truth.index.csv \
    --ref GRCh38_no_alt.fa \
    --outDir ./report \
    --evalType aardvark \
    --batchSystem slurm --binariesMode singularity
```

`--evalType` takes any of `aardvark`, `vcfeval`, `truvari`, `happy`, and more than one at a time.

**Nothing has to be staged by hand.** `--queryVcf`, `--truthIndex` rows and `--ref` may all be
URLs, and `--ref` may be gzipped — it is fetched, decompressed, `samtools faidx`'d and (for
vcfeval/hap.py) turned into an rtg SDF in a single job on a worker. So a whole run can be launched
without downloading a 3 GB fasta first:

```
toil-vcfcomp ./jobstore \
    --queryVcf https://.../hprc-v2.0-mc-chm13.wave.vcf.gz \
    --truthIndex giab.chm13.csv \
    --ref https://hgdownload.soe.ucsc.edu/goldenPath/hs1/bigZips/hs1.fa.gz \
    --outDir ./report --evalType aardvark \
    --batchSystem slurm --binariesMode singularity
```

The reference is still checked against the query and truth VCF headers before any comparison runs:
contig names and lengths must agree, and a mismatch is fatal rather than a warning, because
`bcftools norm -c w` would otherwise pass the mismatched records through and fill the report with
plausible wrong numbers.

### Comparing several call sets

`--queryVcf` takes more than one VCF. Each becomes a **callset**, named after its file, benchmarked
against the same truth and reported side by side — so two graphs can be compared directly instead
of by diffing two separate runs:

```
toil-vcfcomp ./jobstore \
    --queryVcf mc-grch38.vcf.gz mc-grch38.wave.vcf.gz \
    --truthIndex truth.index.csv --ref GRCh38_no_alt.fa --outDir ./report
```

```
callset         sample   tool      variant_type  metric_recall  metric_precision
mc-grch38       HG00097  aardvark  ALL           0.9814         0.9792
mc-grch38.wave  HG00097  aardvark  ALL           0.9841         0.9803
```

`report.summary.tsv` groups by callset as well, so each gets its own across-sample mean/median.
The fan-out is callset × sample × tool, but each sample's truth is fetched **once** and shared by
every callset. Call sets need not have identical sample sets; any sample missing from one of them
is named in the log. Files with the same basename in different directories are disambiguated by
their parent directory.

What it does:

- Splits the query into single-sample VCFs in **one pass** with `bcftools +split`, rather than a
  `bcftools view` per sample over the whole file.
- Builds the rtg SDF **once** when `vcfeval` is requested, instead of once per comparison.
- Works out each sample's sex from its truth VCF -- chrY records are present if and only if the
  sample is male -- and sets the haploid-X / exclude-Y handling (and hap.py's `--gender`)
  accordingly. No sample metadata table needed.
- Restricts each comparison to that sample's own dipcall BED, intersected and merged with
  `--regions` if given.
- Preprocesses **both** sides. `vcfcomp` alone assumes the truth is a ready-to-use GIAB benchmark;
  dipcall output is not, so it gets the same normalization and ploidy handling as the query.
- Applies an allele-length filter to both sides: `--maxLength` (default 1000) for the
  small-variant tools (hap.py, vcfeval) and `--svMaxLength` (default 100000) for the SV-aware ones
  (aardvark, truvari), with `--minLength` (default 0) as the floor. These bound the truth set as
  well as the query, so raise them if you care about alleles longer than that.
- Leaves chrX alone inside the pseudoautosomal regions, where a male is genuinely diploid.
  `--parBed` overrides the built-in GRCh38 PAR coordinates (needed for CHM13).

A sample whose comparison fails does not take the report with it: the tables are written from
whatever succeeded, the failures are listed, and the run is only called a failure if more than
`--maxFailedSamples` (default 0) of them fail.

### Output

```
report/
  report.txt          READ THIS FIRST -- headline numbers, F1 by variant type, and which
                      callset won. a few dozen lines whatever the run size
  report.compare.tsv  the same head-to-head, machine readable: callsets as columns,
                      signed delta against the first --queryVcf, and the winner
  report.tsv          one row per (callset, truthset, sample, tool, comparison, region, filter, variant type)
  report.summary.tsv  the same collapsed across samples, grouped by callset and truthset:
                      n, mean/median/min/max of recall/precision/f1
  report.chrom.tsv    per-contig FP/FN counts, long form
  <callset>/<truthset>/<sample>/<tool>/    that comparison's summary files
```

`report.tsv` is deliberately verbose — aardvark alone emits 45 rows per comparison (five
comparison modes × nine variant types), so a couple of hundred samples is a very large table.
`report.txt` is the digest:

```
VERDICT: hprc-v2.1-mc-chm13.wave has the better F1 on 7 of 10 (tool, truthset, variant type) combinations.

Headline  (variant_type=ALL, mean across samples)
  tool      truthset     callset                  n  recall   precis.  F1
  aardvark  v5.0q_smvar  hprc-v2.0-mc-chm13.wave  1  0.98375  0.98276  0.98326
  aardvark  v5.0q_smvar  hprc-v2.1-mc-chm13.wave  1  0.98502  0.98440  0.98471

F1 by variant type (delta = last callset - baseline)
  truthset     variant_type  hprc-v2.0-...  hprc-v2.1-...  delta      best
  v5.0q_smvar  ALL           0.98326        0.98471        +0.001450  hprc-v2.1-...
  v5.0q_smvar  Insertion     0.97271        0.96982        -0.002890  hprc-v2.0-...
```

It shows each tool's genotype-aware comparison only (`GT` for aardvark/vcfeval/hap.py, `ALLELE`
for truvari, which does not require a genotype match), unstratified by region or filter, so the
numbers in it are like-for-like. Recall and precision sit next to F1 in the headline block so a
precision/recall swap is visible rather than buried inside the F1.

Each tool reports its headline numbers in its own format (aardvark `summary.tsv`, vcfeval
`summary.txt`, truvari `summary.json`, hap.py `*.summary.csv`); all four are normalised onto
aardvark's column set. `--exportAll` additionally exports each comparison's full output directory
as a tarball -- off by default, since the fp/fn/tp VCFs for a few hundred samples add up.

## Worked examples

Both of these compare the same two HPRC minigraph-cactus releases against each other. Nothing is
staged first — every input is a URL, and the reference is fetched, unpacked and indexed on a
worker.

### Against GIAB, on CHM13

GIAB publishes a CHM13v2.0 benchmark for **HG002 only** — the Q100 draft, as separate small-variant
and structural truth sets. So this is 2 call sets × 2 truth sets × 1 sample.

```
giab-truthset --ref CHM13v2.0 --benchmarks all -o giab.chm13.csv

toil-vcfcomp ./jobstore \
    --queryVcf https://s3-us-west-2.amazonaws.com/human-pangenomics/pangenomes/freeze/release2/minigraph-cactus/v2.0/hprc-v2.0-mc-chm13/hprc-v2.0-mc-chm13.wave.vcf.gz \
               https://s3-us-west-2.amazonaws.com/human-pangenomics/pangenomes/freeze/release2/minigraph-cactus/v2.1/hprc-v2.1-mc-chm13/hprc-v2.1-mc-chm13.wave.vcf.gz \
    --truthIndex giab.chm13.csv \
    --ref https://hgdownload.soe.ucsc.edu/goldenPath/hs1/bigZips/hs1.fa.gz \
    --outDir ./giab-chm13-report \
    --evalType aardvark \
    --batchSystem slurm \
    --slurmTime 100:00:00 \
    --doubleMem true \
    --binariesMode singularity \
    --threads 8 --compareMemory 32Gi \
    --logFile toil-vcfcomp.log \
    --logInfo
```

Expect a loud warning that 231 samples are query-only. That is correct here — only HG002 is
shared — but it is the same warning a real sample-name mismatch produces, so do not learn to
ignore it.

### Against the HPRC dipcall calls, on CHM13

The Yale R2 submission has dipcall calls on CHM13v2 as well as GRCh38, which gets you **231
samples** instead of one. Same graphs, same reference, only the manifest changes:

```
hprc-truthset --ref CHM13v2 --caller dipcall -o hprc.chm13.csv

toil-vcfcomp ./jobstore-hprc \
    --queryVcf https://s3-us-west-2.amazonaws.com/human-pangenomics/pangenomes/freeze/release2/minigraph-cactus/v2.0/hprc-v2.0-mc-chm13/hprc-v2.0-mc-chm13.wave.vcf.gz \
               https://s3-us-west-2.amazonaws.com/human-pangenomics/pangenomes/freeze/release2/minigraph-cactus/v2.1/hprc-v2.1-mc-chm13/hprc-v2.1-mc-chm13.wave.vcf.gz \
    --truthIndex hprc.chm13.csv \
    --ref https://hgdownload.soe.ucsc.edu/goldenPath/hs1/bigZips/hs1.fa.gz \
    --outDir ./hprc-chm13-report \
    --evalType aardvark \
    --batchSystem slurm \
    --slurmTime 100:00:00 \
    --doubleMem true \
    --binariesMode singularity \
    --threads 8 --compareMemory 32Gi \
    --maxJobs 200 \
    --logFile toil-vcfcomp-hprc.log \
    --logInfo
```

That is **462 whole-genome comparisons** (2 call sets × 231 samples), so unlike the GIAB run it is
a real cluster job. `--maxJobs` keeps it from taking the whole queue. Each sample's truth is
downloaded once and shared by both call sets, so it is 231 downloads rather than 462.

A note on the graphs' own samples: each carries 232, of which 231 have dipcall truth. The
remaining one is `GRCh38`, a reference path in the graph rather than a sample, and it is correctly
reported as query-only.

## Things that will bite you

**`GT="alt"` is not the filter you want.** bcftools classifies a half-called genotype like `.|1`
as `"mis"`, not `"alt"`, so filtering a graph VCF with `GT="alt"` silently drops every hemizygous
call: all of chrY and non-PAR chrX in males, and around 40% of calls overall in a typical
deconstructed graph VCF where one haplotype often does not traverse a site. `--splitFilter`
defaults to `GT~"[1-9]"`, which keeps any genotype carrying a non-ref allele.

**dipcall's good calls are `FILTER="."`, not `PASS`.** Only about 4% of records carry a filter
(`HET1`/`HET2`/`GAP1`/`GAP2`/`DIPX`). A `--pass-only` style filter therefore discards roughly 96%
of the truth set and reports a recall near zero. `--truthFilter` defaults to
`FILTER="." || FILTER="PASS"`. hap.py's stock options include `--pass-only`, so the workflow warns
if you ask for `happy` without removing it.

**`--workDir` has to be somewhere the container runtime can read.** A snap-confined docker
bind-mounts an *empty* directory for any host path outside `$HOME`, with no error -- the first
symptom is a tool reporting that its inputs do not exist. The workflow probes for this up front
and tells you to move `--workDir`; Toil otherwise defaults to `/tmp`.

**Singularity on a cluster.** `--binariesMode singularity` converts each docker image into a
sandbox directory under `$SINGULARITY_CACHEDIR/toil-vcfcomp`, built once and renamed into place
atomically, because concurrent `singularity pull`s sharing a cache corrupt each other. The images
are prepared in a single job before the fan-out. The work directory is bind-mounted at `/mnt`,
not `/data`: only one of the four tool images ships a `/data`, and a singularity install without
overlay support silently skips a bind whose destination does not exist.

When `--binariesMode` is not given it is autodetected, preferring singularity whenever
`--batchSystem` is not `single_machine` — the submit node commonly has docker while the compute
nodes do not, so probing for docker there would pick a runtime no worker can use.

**The reference must match.** `bcftools norm -c w` passes a reference-mismatched record through
unchanged and exits 0, so the wrong `--ref` yields a completed run full of plausible, wrong
numbers. The reference job compares contig names and lengths against the query and truth VCF
headers and refuses to continue on a mismatch, before any comparison runs.

## Sizing

Per-comparison resources are `--threads` cores, `--compareMemory` (default 32Gi) and
`5 x <reference size> + --compareDisk` (default 20Gi) of disk. Slurm enforces memory through
cgroups, so lower these if your jobs are being killed and raise them if aardvark runs out.
