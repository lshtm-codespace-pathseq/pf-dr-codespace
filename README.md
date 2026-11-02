<!--
  To enable the one-click launch button, replace OWNER/REPO below with your
  GitHub repo (e.g. janedoe/pf-dr-workshop), then uncomment the line.
  [![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/OWNER/REPO)
-->

# Detecting malaria drug resistance from genome sequence — a 2-hour hands-on course

Welcome! Over the next two hours you will take **raw DNA sequencing reads from a
single *Plasmodium falciparum* (malaria) sample** and work, one command at a time,
all the way to a **drug-resistance report** — the kind of result that tells a
public-health team which antimalarials a parasite population is likely to resist.

You will run every step yourself and see exactly what each tool does:

```
FastQC  →  trim (fastp)  →  remove human reads (hostile)  →  map to genome (bwa mem -M)
       →  mark duplicates  →  coverage (mosdepth)  →  call variants (GATK)
       →  hard-filter  →  interpret resistance markers  →  plot mdr1 copy number  →  MultiQC
```

Everything runs **in your web browser** through GitHub Codespaces — there is
**nothing to install** on your own computer. The tools are pre-installed for you.

> **This is a teaching pipeline, for training and research — not for clinical
> diagnosis.** For real surveillance, see the curated resources linked at the end.

---

## Learning objectives

By the end of this course you will be able to:

1. **Explain the whole short-read variant-calling workflow** and why each stage
   exists — QC, trimming, host removal, mapping, duplicate marking, coverage,
   variant calling, filtering.
2. **Run each step as a real command** and read its output to judge whether it
   worked (quality scores, mapping rate, duplicate rate, depth, filter status).
3. **Interpret the biology**: translate variants into amino-acid changes at known
   resistance codons, and recognise **gene amplification** (copy number) as a
   second, separate resistance mechanism.
4. **Name the main *P. falciparum* resistance genes** (*crt, mdr1, dhfr, dhps,
   kelch13, cytb*) and the drugs they affect.
5. **Describe how this scales** from one sample to a cohort, and when to swap
   hard-filtering for machine-learning filtering (VQSR).

**Who it's for:** anyone comfortable typing commands into a terminal. No prior
bioinformatics needed — every command is given and explained. Some familiarity
with DNA, genes and the idea of a "mutation" helps.

---

## Background: why this matters (2-minute read)

Malaria killed an estimated 597,000 people in 2023, most of them children in
sub-Saharan Africa ([WHO World Malaria Report](https://www.who.int/teams/global-malaria-programme/reports/world-malaria-report-2024)).
*P. falciparum* is the deadliest species. Our main weapons are **drugs** — and the
parasite keeps evolving **resistance** to them. Chloroquine, the former front-line
drug, was lost to resistance across most of the world; **artemisinin** (the core of
today's first-line ACT combinations) is now showing **partial resistance** in
Southeast Asia and, more recently, East Africa.

Resistance is driven by **changes in a handful of parasite genes**. Two kinds of
change matter:

- **Point mutations** (single-letter DNA changes) that alter a protein — e.g.
  *crt* **K76T** for chloroquine, *kelch13* **C580Y** for artemisinin.
- **Gene amplification** (extra *copies* of a gene) — e.g. extra copies of **mdr1**
  reduce susceptibility to lumefantrine and mefloquine.

If we **sequence** a parasite's genome, we can read these changes directly. That
lets surveillance programmes track resistance across a country and choose drugs
that still work. This course is a miniature version of exactly that analysis.

Good primers if you want more depth:
- Review: *Antimalarial drug resistance: linking P. falciparum parasite biology to
  the clinic* — Blasco, Leroy & Fidock, **Nature Medicine** 2017,
  [doi:10.1038/nm.4381](https://doi.org/10.1038/nm.4381).
- [WHO Malaria Threats Map](https://apps.who.int/malaria/maps/threats/) — live view
  of where resistance markers are found.
- [MalariaGEN](https://www.malariagen.net/) — the global *P. falciparum* genomic
  epidemiology project whose methods this course mirrors.

---

## Our example sample

**ERR018972** — a paired-end Illumina whole-genome sequencing run of *P.
falciparum*, hosted on the **European Nucleotide Archive (ENA)**. You can see its
record here: <https://www.ebi.ac.uk/ena/browser/view/ERR018972>.

We map it to the **3D7 reference genome** — the standard lab reference strain —
downloaded from **[PlasmoDB](https://plasmodb.org/)**, the *Plasmodium* genomics
database. The 3D7 genome is ~23 Mb across 14 chromosomes plus the mitochondrion and
apicoplast, and is famously **AT-rich (~81% A+T, only ~19% GC)** — remember that
when the GC plots look unusual; it's the biology, not a problem.

---

## Part 1 — Start your Codespace (≈ 5 min)

A **Codespace** is a full Linux computer that runs in the cloud and opens as a
VS Code editor **inside your browser**. GitHub builds it from this repository and
installs all the tools automatically. Nothing lands on your laptop.

1. On this repository's page on GitHub, click the green **`< > Code`** button ▸
   **Codespaces** tab ▸ **Create codespace on main**.
2. A browser editor opens. **The first time, it spends ~3–5 minutes building the
   environment** (installing FastQC, fastp, hostile, bwa, samtools, GATK, mosdepth,
   MultiQC …). You'll see progress in the terminal; it finishes with
   `Setup complete`. *(If your instructor enabled "prebuilds", this is near-instant.)*
3. **Open a fresh terminal**: menu ▸ **Terminal ▸ New Terminal**. The prompt should
   start with **`(pfdr)`** — that means the software environment is active.

**Check you're ready** — both of these should print something sensible:

```bash
gatk --version
bwa 2>&1 | head -3
```

<details>
<summary>Prompt doesn't show <code>(pfdr)</code>?</summary>

The setup may still be running (wait for `Setup complete`), or just activate it:

```bash
eval "$("$HOME/.local/bin/micromamba" shell hook -s bash)" && micromamba activate pfdr
```
</details>

> **New to the terminal?** You type a command and press Enter; it runs and prints
> output. `mkdir -p X` makes a folder, `ls X` lists it, `|` sends one tool's output
> into the next. That's 90% of what you need today.

---

## Part 2 — Run the pipeline, one command at a time

Work through the steps **in order**. Read the output of each before moving on —
noticing *what changed* is where the learning is.

First, set two variables you'll reuse everywhere (a Codespace has ~4 CPU cores):

```bash
S=ERR018972     # the sample name
T=4             # number of CPU threads to use
```

> **Steps 00 and 01 are just downloads**, so they stay as one-line scripts. **Every
> other step you type out as the real commands** — that's the point of the course.
> The identical commands also live in `steps/NN_*.sh` if you want to check your
> typing, and `bash run_all.sh` runs the whole thing start to finish.
>
> **Short on time?** See *Workshop mode* right after step 01 to make everything run
> faster.

---

### Step 00 — Get the reference genome *(script)*

```bash
mkdir ref
cd ref
wget -O Pf3D7.fasta "https://zenodo.org/records/18681980/files/PlasmoDB-54-Pfalciparum3D7-Genome.fasta?download=1"
wget -O Pf3D7.gff.gz "https://zenodo.org/records/18681980/files/PlasmoDB-55_Pfalciparum3D7.gff.gz?download=1"
gunzip Pf3D7.gff

# Build indexes
samtools faidx Pf3D7.fasta
bwa index Pf3D7.fasta
gatk CreateSequenceDictionary -R Pf3D7.fasta

python ./scripts/build_dr_regions.py' \
      --fasta Pf3D7.fasta \
      --gff Pf3D7.fasta.gff \
      --genes ./resources/dr_genes.tsv \
      --out-bed ./resources/dr_genes.bed \
      --out-cds ./resources/dr_genes.cds.json

cd ../
```

**What it does.** Mapping and variant-calling tools can't just read a plain FASTA —
they need pre-built **index** files beside it. This script downloads the 3D7 genome
from PlasmoDB and builds:

- a **`.fai`** index (`samtools faidx`) — lets tools jump to any position instantly;
- a **`.dict`** dictionary (`gatk CreateSequenceDictionary`) — lists the contigs;
- the **BWA index** (`bwa index`) — the data structure the aligner searches against;
- from the PlasmoDB **GFF** (gene annotation), a **BED** of the drug-resistance
  genes and a **CDS model** (`ref/dr_genes_cds.json`) used later to turn DNA
  variants into amino-acid changes.

**Why pre-indexing?** A genome is large; searching it read-by-read without an index
would be hopelessly slow. The index is built once and reused by every later step —
the same idea as the index at the back of a textbook.

> **Check:** the command ends by listing files in `ref/`. You should see
> `Pf3D7.fasta`, `Pf3D7.fasta.bwt` (and friends), `Pf3D7.dict`, and `dr_genes.bed`.

### Step 01 — Download the sequencing reads *(script)*

```bash
bash steps/01_get_reads.sh
```

**What it does.** Downloads the two FASTQ files for our sample from the ENA:
`data/${S}_1.fastq.gz` (forward reads, "R1") and `data/${S}_2.fastq.gz` (reverse
reads, "R2"). **Paired-end** sequencing reads each DNA fragment from both ends,
which helps place reads accurately on the genome.

<details>
<summary>What's a FASTQ file?</summary>

Four lines per read: an ID, the DNA sequence, a `+`, and a **quality string** —
one character per base encoding how confident the sequencer was (the "Phred"
score). FastQC (next step) summarises those qualities. Spec and background:
[Illumina FASTQ](https://support.illumina.com/bulletins/2016/04/fastq-files-explained.html).
</details>

#### ⏩ Workshop mode (recommended for a timed class)

The full read set takes a while to download and process. To keep the session
snappy, **subsample** to e.g. 2 million read pairs:

```bash
SUBSAMPLE_PAIRS=2000000 bash steps/01_get_reads.sh
```

Every later step auto-detects the subsampled files — **no other change needed** if
you use the `steps/*.sh` scripts. If you're typing the raw commands below, swap
`data/${S}_1.fastq.gz` → `data/${S}.sub_1.fastq.gz` (and `_2`) in steps 02–03.

> Trade-off: fewer reads = **lower coverage**, which makes the mdr1 copy-number plot
> (step 11) and mixed-infection calls less clear. Fine for learning; use the full
> set for a "real" result.

---

### Step 02 — Look at the raw data with FastQC

**Golden rule of sequencing analysis: look at your data before you touch it.**
[FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) reads a FASTQ
and produces an HTML report of quality metrics.

```bash
mkdir -p results/01_fastqc_raw

fastqc -t $T --outdir results/01_fastqc_raw \
    data/${S}_1.fastq.gz data/${S}_2.fastq.gz
```

**Open** `results/01_fastqc_raw/${S}_1_fastqc.html` (right-click in the file
explorer ▸ *Open with Live Preview*, or download it).

**What to look at:**
- **Per-base sequence quality** — a box plot of quality along the read. It usually
  dips toward the 3′ end; that's normal and we'll trim it next.
- **Adapter content** — leftover sequencing adapter, also trimmed next.
- **Per-sequence GC content** — for pure *P. falciparum* this peaks very low
  (~19% GC). A **second peak around ~41%** is a tell-tale sign of **human DNA
  contamination** — which step 04 removes.

Guide to each FastQC module (what "good" vs "bad" looks like):
[FastQC Help](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/).

> **✅ Check your understanding**
> 1. Why do quality scores typically fall off at the end of a read?
> 2. If you saw a strong GC peak at ~41%, what would you suspect — and which later
>    step addresses it?

---

### Step 03 — Trim adapters and low-quality bases (fastp), then re-check

[fastp](https://github.com/OpenGene/fastp) cleans reads: it removes sequencing
**adapters**, trims **low-quality** bases from the ends, and discards reads that
become too short. `--detect_adapter_for_pe` auto-detects adapters for paired-end
data. We then run FastQC again so you can compare **before vs after**.

```bash
mkdir -p results/02_trim

fastp \
    --in1 data/${S}_1.fastq.gz --in2 data/${S}_2.fastq.gz \
    --out1 results/02_trim/${S}.trim_1.fastq.gz \
    --out2 results/02_trim/${S}.trim_2.fastq.gz \
    --detect_adapter_for_pe \
    --cut_front --cut_tail --cut_mean_quality 20 \
    --qualified_quality_phred 20 \
    --length_required 50 \
    --thread $T \
    --json results/02_trim/fastp.json \
    --html results/02_trim/fastp.html

fastqc -t $T --outdir results/02_trim \
    results/02_trim/${S}.trim_1.fastq.gz \
    results/02_trim/${S}.trim_2.fastq.gz
```

**What the flags mean:** `--cut_front --cut_tail --cut_mean_quality 20` slides along
each read and trims ends whose average quality drops below Q20 (= 99% base-call
accuracy); `--qualified_quality_phred 20` defines what counts as a "good" base;
`--length_required 50` drops reads shorter than 50 bp after trimming.

**Open** `results/02_trim/fastp.html` — a single, clear before/after summary.

> **✅ Check your understanding**
> - A Phred quality of 20 means a 1-in-100 chance the base is wrong. What error
>   probability does **Q30** mean? *(Hint: Phred Q = −10·log₁₀(error).)*
>
> **🧪 Exercise:** re-run fastp with a stricter `--cut_mean_quality 30` into a
> different output folder and compare how many reads survive in `fastp.html`.

---

### Step 04 — Remove human (host) reads with hostile

Clinical malaria samples are sequenced from **patient blood**, so the DNA is a
**mixture of parasite and human**. We must remove human reads before analysis:
it speeds things up, avoids false alignments, and — importantly — means no human
genetic data is left in files we might share.
[hostile](https://github.com/bede/hostile) does this by aligning reads to a human
reference and **keeping only the reads that do *not* match**.

```bash
mkdir -p results/03_dehost

hostile clean \
    --fastq1 results/02_trim/${S}.trim_1.fastq.gz \
    --fastq2 results/02_trim/${S}.trim_2.fastq.gz \
    --aligner bowtie2 \
    --index human-t2t-hla \
    --threads $T \
    --out-dir results/03_dehost \
    --force > results/03_dehost/hostile_log.json

# hostile names its outputs *.clean_1.fastq.gz / *.clean_2.fastq.gz —
# copy them to tidy names the next step expects:
cp results/03_dehost/*clean_1.fastq.gz results/03_dehost/${S}.clean_1.fastq.gz
cp results/03_dehost/*clean_2.fastq.gz results/03_dehost/${S}.clean_2.fastq.gz
```

> ⏳ **First run downloads a ~3.5 GB human genome index** (`human-t2t-hla`), then
> caches it — so this step is the slow one. Start it, then read on while it works.
> (Your instructor may have pre-cached it, making this instant.)

**What to look at:** `results/03_dehost/hostile_log.json` reports `reads_in`,
`reads_removed` and `reads_removed_proportion` — the fraction flagged as human.
hostile deliberately errs toward **keeping** microbial reads. Method and paper:
Constantinides *et al.*, **Bioinformatics** 2023,
[doi:10.1093/bioinformatics/btad728](https://doi.org/10.1093/bioinformatics/btad728).

> **✅ Check your understanding:** why is it scientifically *and* ethically
> important to strip human reads from a pathogen sequencing dataset?

---

### Step 05 — Map the reads to the genome with `bwa mem -M`

Now we find **where each read came from** on the 3D7 genome.
[bwa mem](https://github.com/lh3/bwa) is the standard short-read aligner.

```bash
mkdir -p results/04_map

bwa mem -M -t $T \
    -R "@RG\tID:${S}\tSM:${S}\tLB:${S}_lib1\tPL:ILLUMINA" \
    ref/Pf3D7.fasta \
    results/03_dehost/${S}.clean_1.fastq.gz \
    results/03_dehost/${S}.clean_2.fastq.gz \
  | samtools sort -@ $T -o results/04_map/${S}.sorted.bam -

samtools index results/04_map/${S}.sorted.bam
samtools flagstat results/04_map/${S}.sorted.bam
```

**The two flags that matter:**
- **`-M`** marks shorter "split" alignments as *secondary*. This keeps the output
  compatible with Picard/GATK, which expect exactly one primary alignment per read.
  *(This is the flag the course specifies.)*
- **`-R`** attaches a **read group** header (`@RG`). It records which **sample**
  (`SM`) and **library** (`LB`) the reads came from. **GATK refuses to run without
  it** — on a real project the read group is how reads from many lanes/samples stay
  correctly labelled. Background:
  [GATK read groups](https://gatk.broadinstitute.org/hc/en-us/articles/360035890671-Read-groups).

We then pipe (`|`) bwa's output straight into **`samtools sort`** (no giant
intermediate file) to make a coordinate-sorted **BAM** (the compressed binary form
of an alignment), and **index** it. [samtools](https://www.htslib.org/) is the
Swiss-army knife for BAM files.

**What to look at:** the **`flagstat`** summary — especially the **mapped %** (what
fraction of reads found a home on the genome) and **properly paired %**. Decode any
SAM flag here: [Explain SAM flags](https://broadinstitute.github.io/picard/explain-flags.html).

> **🧪 Exercise:** peek at the actual alignments —
> `samtools view results/04_map/${S}.sorted.bam | head -3` — and find the mapping
> position (column 4) and mapping quality (column 5).

---

### Step 06 — Mark duplicate reads

During library prep, PCR can copy the **same original DNA fragment** many times.
Those duplicate read-pairs aren't independent evidence, so we must not let each one
"vote" during variant calling. `MarkDuplicates` finds them and **flags** them in
the BAM (it doesn't delete them); the caller then ignores flagged reads.

```bash
mkdir -p results/05_markdup

gatk MarkDuplicates \
    -I results/04_map/${S}.sorted.bam \
    -O results/05_markdup/${S}.markdup.bam \
    -M results/05_markdup/${S}.dup_metrics.txt \
    --CREATE_INDEX true
```

**What to look at:** `PERCENT_DUPLICATION` in `${S}.dup_metrics.txt`. High
duplication means low library complexity (little unique DNA). Tool docs:
[MarkDuplicates](https://gatk.broadinstitute.org/hc/en-us/articles/360037052812-MarkDuplicates-Picard).

> **A note on BQSR.** A human-genome pipeline would next run **Base Quality Score
> Recalibration** (`BaseRecalibrator` + `ApplyBQSR`) to correct systematic errors
> in the quality scores. BQSR needs a trusted set of *known variants* to mask real
> variation — and *P. falciparum* has no single standard known-sites set, so this
> course **skips BQSR**. See [Scaling to a cohort](#scaling-to-a-cohort) for how
> you'd add it. Docs:
> [BQSR](https://gatk.broadinstitute.org/hc/en-us/articles/360035890531-Base-Quality-Score-Recalibration-BQSR).

---

### Step 07 — Measure coverage with mosdepth

**Coverage** (or "depth") = how many reads overlap each position. We need it for two
reasons: to judge whether we have *enough* data to trust a variant, and to spot
**gene amplification** (a duplicated gene has roughly double the surrounding depth).
[mosdepth](https://github.com/brentp/mosdepth) computes it fast.

```bash
mkdir -p results/06_coverage

# (a) genome-wide, in 500 bp windows — the track we'll plot across mdr1
mosdepth --no-per-base --fast-mode -t $T --by 500 \
    results/06_coverage/${S}.genome \
    results/05_markdup/${S}.markdup.bam

# (b) restricted to the drug-resistance genes — a quick per-gene table
mosdepth --no-per-base -t $T --by ref/dr_genes.bed \
    results/06_coverage/${S}.dr \
    results/05_markdup/${S}.markdup.bam

# peek at the per-gene mean depth:
zcat results/06_coverage/${S}.dr.regions.bed.gz
```

**What to look at:** in the per-gene table, is **mdr1** noticeably deeper than the
other genes? That's the amplification signal you'll visualise in step 11. Check the
genome-wide average with:

```bash
awk '$1=="total"' results/06_coverage/${S}.genome.mosdepth.summary.txt
```

mosdepth paper: Pedersen & Quinlan, **Bioinformatics** 2018,
[doi:10.1093/bioinformatics/btx699](https://doi.org/10.1093/bioinformatics/btx699).

> **✅ Check your understanding:** why is coverage, not just the presence of a
> variant read, essential before you believe a genotype call?

---

### Step 08 — Call variants with GATK (the GVCF way)

Now we find the positions where our sample **differs from the reference**.
We use [GATK](https://gatk.broadinstitute.org/)'s **HaplotypeCaller**, which
locally re-assembles reads to call SNPs and short indels accurately. This is done
in **two moves**, exactly as in the
[GATK germline best-practice workflow](https://gatk.broadinstitute.org/hc/en-us/articles/360035535932-Germline-short-variant-discovery-SNPs-Indels):

```bash
mkdir -p results/07_variants

# 1) HaplotypeCaller in GVCF mode: a record for EVERY position
gatk HaplotypeCaller \
    -R ref/Pf3D7.fasta \
    -I results/05_markdup/${S}.markdup.bam \
    -O results/07_variants/${S}.g.vcf.gz \
    -ERC GVCF \
    --sample-ploidy 2 \
    -L ref/dr_genes.bed --interval-padding 1000

# 2) GenotypeGVCFs: turn the GVCF into a final genotyped VCF
gatk GenotypeGVCFs \
    -R ref/Pf3D7.fasta \
    -V results/07_variants/${S}.g.vcf.gz \
    -O results/07_variants/${S}.vcf.gz

# how many variants were called?
zcat results/07_variants/${S}.vcf.gz | grep -vc '^#'
```

**Why two moves for one sample?** The intermediate **GVCF**
([what's a GVCF?](https://gatk.broadinstitute.org/hc/en-us/articles/360035531812-GVCF-Genomic-Variant-Call-Format))
records *every* position, variant or not. For one sample it seems like extra work —
but it's exactly what lets you **scale to a cohort**: make one GVCF per sample, then
genotype them all together (see [below](#scaling-to-a-cohort)). We teach the
scalable pattern from the start. HaplotypeCaller docs:
[HaplotypeCaller](https://gatk.broadinstitute.org/hc/en-us/articles/360037225632-HaplotypeCaller).

**`-L ref/dr_genes.bed`** restricts calling to the drug-resistance genes so it
finishes quickly. To call the **whole genome** instead, drop that line (slower).

**Ploidy.** Blood-stage *P. falciparum* is **haploid** (one copy of each
chromosome). But a patient can carry **more than one parasite clone** — a *mixed
infection*. Calling with **`--sample-ploidy 2`** lets a second clone show up as a
"heterozygous" genotype, which is how [MalariaGEN](https://www.malariagen.net/)
handles it. The report in step 10 flags such sites as *mixed*.

> **✅ Check your understanding:** a true, pure haploid sample should have no
> heterozygous calls. If you see many, what might that tell you about the sample?

---

### Step 09 — Hard-filter the variants

Not every called variant is real; some are artefacts of mapping or sequencing. We
flag the suspicious ones using fixed thresholds on quality annotations — **"hard
filtering"**. It's the right tool when you *don't* have enough data to train the
machine-learning alternative (VQSR). SNPs and indels misbehave differently, so we
**split → filter each → merge**. Filtered records are **labelled** in the FILTER
column (not deleted), so you can always see *why* something was dropped.

```bash
mkdir -p results/08_filter

# 1) split into SNPs and indels
gatk SelectVariants -R ref/Pf3D7.fasta -V results/07_variants/${S}.vcf.gz \
    --select-type-to-include SNP   -O results/08_filter/${S}.snps.raw.vcf.gz
gatk SelectVariants -R ref/Pf3D7.fasta -V results/07_variants/${S}.vcf.gz \
    --select-type-to-include INDEL -O results/08_filter/${S}.indels.raw.vcf.gz

# 2a) SNP hard filters (GATK-recommended starting thresholds)
gatk VariantFiltration -R ref/Pf3D7.fasta \
    -V results/08_filter/${S}.snps.raw.vcf.gz \
    -O results/08_filter/${S}.snps.filt.vcf.gz \
    --filter-expression "QD < 2.0"              --filter-name "QD2" \
    --filter-expression "FS > 60.0"             --filter-name "FS60" \
    --filter-expression "MQ < 40.0"             --filter-name "MQ40" \
    --filter-expression "MQRankSum < -12.5"     --filter-name "MQRankSum" \
    --filter-expression "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum" \
    --filter-expression "SOR > 3.0"             --filter-name "SOR3"

# 2b) indel hard filters (looser strand/position thresholds)
gatk VariantFiltration -R ref/Pf3D7.fasta \
    -V results/08_filter/${S}.indels.raw.vcf.gz \
    -O results/08_filter/${S}.indels.filt.vcf.gz \
    --filter-expression "QD < 2.0"               --filter-name "QD2" \
    --filter-expression "FS > 200.0"             --filter-name "FS200" \
    --filter-expression "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum" \
    --filter-expression "SOR > 10.0"             --filter-name "SOR10"

# 3) merge back together
gatk MergeVcfs \
    -I results/08_filter/${S}.snps.filt.vcf.gz \
    -I results/08_filter/${S}.indels.filt.vcf.gz \
    -O results/08_filter/${S}.filtered.vcf.gz

# PASS vs flagged breakdown
zcat results/08_filter/${S}.filtered.vcf.gz | grep -v '^#' | cut -f7 | sort | uniq -c
```

**What to look at:** the final count of `PASS` vs flagged records. The DR report in
step 10 keeps only `PASS`. The annotations (QD, FS, MQ, SOR, …) are explained in the
[glossary below](#annotation-glossary). GATK's own walkthrough:
[Hard-filtering germline short variants](https://gatk.broadinstitute.org/hc/en-us/articles/360035890471-Hard-filtering-germline-short-variants).
The VCF format itself: [VCF spec / hts-specs](https://samtools.github.io/hts-specs/).

> **✅ Check your understanding:** why do indels get *looser* strand-bias and
> read-position thresholds than SNPs? *(Hint: think about how hard indels are to
> align near read ends.)*

---

### Step 10 — Interpret the resistance markers 🎯

This is the payoff. The helper script `scripts/report_dr.py` takes each known
resistance marker (gene, codon, wild-type and resistant amino acids) and:

1. finds that codon's position in the genome from the CDS model (**correctly
   handling genes on the + and − DNA strands**),
2. **self-checks** that the reference codon really encodes the documented
   wild-type amino acid (a guard against a reference/annotation version mismatch),
3. applies the PASS variants, **translates** the codon, and classifies it as
   **wild-type / resistance marker / other mutation**,
4. reports allele fractions so **mixed infections** are visible.

```bash
mkdir -p results/09_report

python scripts/report_dr.py \
    --fasta ref/Pf3D7.fasta \
    --cds ref/dr_genes_cds.json \
    --markers resources/dr_markers.tsv \
    --vcf results/08_filter/${S}.filtered.vcf.gz \
    --coverage results/06_coverage/${S}.dr.regions.bed.gz \
    --sample $S \
    --outdir results/09_report \
    --min-dp 5 \
    --mixed-min-af 0.10
```

It prints a table and writes `results/09_report/${S}.dr_variants.tsv` (every
marker), `${S}.dr_calls.tsv` (only the non-wild-type calls), and a JSON.

**The marker list lives in [`resources/dr_markers.tsv`](resources/dr_markers.tsv)**
— a plain table you can read and extend. It covers *crt, mdr1, dhfr, dhps, kelch13*
and *cytb*. It is a **teaching set, not a clinical-grade panel** (see
[Going further](#going-further) for curated tools).

> **✅ Check your understanding:** your report shows **crt K76T**. Which drug does
> that marker most classically predict resistance to?
>
> **🧪 Exercise:** open `resources/dr_markers.tsv`, add a *kelch13* codon not
> already listed (e.g. **P553L**), and re-run step 10. Did a new row appear?

---

### Step 11 — Plot mdr1 coverage and estimate copy number 📈

*mdr1* resistance has **two** mechanisms: point mutations (handled above) **and
gene amplification**. Amplification can't be seen in the VCF — you see it in the
**depth**. We estimate:

```
copy number  ≈  median depth over mdr1  /  genome-wide median depth
```

```bash
mkdir -p results/10_plots

python scripts/plot_coverage.py \
    --windows results/06_coverage/${S}.genome.regions.bed.gz \
    --summary results/06_coverage/${S}.genome.mosdepth.summary.txt \
    --dr-regions results/06_coverage/${S}.dr.regions.bed.gz \
    --cds ref/dr_genes_cds.json \
    --sample $S \
    --gene mdr1 --gene-id PF3D7_0523000 \
    --out results/10_plots/${S}.mdr1_coverage.png
```

**Open** `results/10_plots/${S}.mdr1_coverage.png`. The top panel is the depth
profile across *mdr1* with the 1× and 2× baselines drawn in; the bottom panel is
mean depth per DR gene. A value clearly **above 1** means amplification, linked to
**reduced lumefantrine and mefloquine susceptibility**. The *same logic* applies to
**plasmepsin 2/3** amplification and **piperaquine** resistance.

> **🧪 Exercise:** copy the command but point it at plasmepsin 2
> (`--gene plasmepsin2 --gene-id PF3D7_1408000`) to plot a different gene's depth.

---

### Step 12 — Pull it all together with MultiQC

[MultiQC](https://multiqc.info/) scans the whole `results/` folder and stitches
every report (FastQC, fastp, samtools, mosdepth, GATK) into **one interactive HTML
page** — the standard way to review a run, and invaluable when comparing many
samples.

```bash
mkdir -p results/11_multiqc

multiqc -f -o results/11_multiqc -n ${S}_multiqc.html results
```

**Open** `results/11_multiqc/${S}_multiqc.html` — your whole analysis on one page.
Paper: Ewels *et al.*, **Bioinformatics** 2016,
[doi:10.1093/bioinformatics/btw354](https://doi.org/10.1093/bioinformatics/btw354).

🎉 **You've gone from raw reads to a drug-resistance report.**

---

## What a finished run looks like

After all steps, `results/` contains:

```
01_fastqc_raw/   raw-read FastQC reports
02_trim/         fastp output + trimmed-read FastQC + fastp.html
03_dehost/       human-free reads + hostile_log.json
04_map/          sorted, indexed BAM + flagstat
05_markdup/      duplicate-marked BAM + metrics
06_coverage/     mosdepth genome-wide + per-DR-gene depth
07_variants/     GVCF + genotyped VCF
08_filter/       hard-filtered VCF (PASS/flagged)
09_report/       dr_variants.tsv, dr_calls.tsv, dr_report.json   ← the answer
10_plots/        <sample>.mdr1_coverage.png                      ← the plot
11_multiqc/      <sample>_multiqc.html                           ← the overview
```

The headline outputs are **`09_report/${S}.dr_variants.tsv`** and
**`10_plots/${S}.mdr1_coverage.png`**.

---

## The drug-resistance biology in one table

| Gene | PlasmoDB ID | Mechanism | Drugs affected |
|------|-------------|-----------|----------------|
| **crt** | PF3D7_0709000 | point mutations (K76T …; plus piperaquine-associated SNPs) | chloroquine, amodiaquine, piperaquine |
| **mdr1** | PF3D7_0523000 | point mutations **and copy number** | lumefantrine, mefloquine, amodiaquine, chloroquine |
| **dhfr** | PF3D7_0417200 | point mutations (N51I, C59R, S108N, I164L) | pyrimethamine (part of SP) |
| **dhps** | PF3D7_0810800 | point mutations (A437G, K540E, A581G) | sulfadoxine (part of SP) |
| **kelch13** | PF3D7_1343700 | propeller-domain mutations (e.g. C580Y, R561H) | artemisinin (partial resistance) |
| **cytb** | PF3D7_MIT02300 | mutation at codon 268 | atovaquone |
| **plasmepsin 2/3** | PF3D7_1408000 / …8100 | copy number | piperaquine |

Full marker-by-marker detail: [`resources/dr_markers.tsv`](resources/dr_markers.tsv).
Foundational references: *crt* K76T as a chloroquine-resistance marker — Djimdé
*et al.*, **NEJM** 2001, [doi:10.1056/NEJM200101253440403](https://doi.org/10.1056/NEJM200101253440403);
*kelch13* and artemisinin resistance — Ariey *et al.*, **Nature** 2014,
[doi:10.1038/nature12876](https://doi.org/10.1038/nature12876). Current validated
*kelch13* markers are maintained by the [WHO](https://www.who.int/teams/global-malaria-programme).

---

## Scaling to a cohort

This course does **one** sample. A real surveillance study does hundreds. The GVCF
design (step 08) is what makes that work: one GVCF per sample, then **joint
genotyping** so every position is evaluated in every sample.

```bash
# 1) per sample (in parallel): steps 00–08 -> sample1.g.vcf.gz, sample2.g.vcf.gz, ...

# 2) consolidate the GVCFs
gatk GenomicsDBImport \
     -V sample1.g.vcf.gz -V sample2.g.vcf.gz ... \
     --genomicsdb-workspace-path cohort_db \
     -L ref/dr_genes.bed

# 3) joint genotyping across the whole cohort
gatk GenotypeGVCFs -R ref/Pf3D7.fasta -V gendb://cohort_db -O cohort.vcf.gz
```

### Why this course skips BQSR and VQSR

- **BQSR** (base quality recalibration) needs a trusted *known-variants* set to
  mask real variation; there isn't a single standard one for *P. falciparum*, and
  the benefit for teaching is small. For production you'd derive a high-confidence
  site list (e.g. from [MalariaGEN Pf7](https://www.malariagen.net/parasite/pf7))
  and run `BaseRecalibrator`/`ApplyBQSR` between steps 06 and 08.
- **VQSR** (machine-learning variant filtering) *replaces* step 09's hard filters,
  but it must be **trained on many samples** with good truth sets — a **cohort**
  tool, not a single-sample one. With a cohort and training resources, swap step 09
  for `VariantRecalibrator` + `ApplyVQSR`. Until then, hard-filtering is correct.
  Docs: [VQSR](https://gatk.broadinstitute.org/hc/en-us/articles/360035531612-Variant-Quality-Score-Recalibration-VQSR).

---

## Annotation glossary

The quality annotations used in step 09's hard filters (GATK-recommended starting
thresholds):

| Annotation | Plain meaning | SNP filter | Indel filter |
|-----------|---------------|-----------|--------------|
| **QD** | variant quality divided by depth (quality per read) | `< 2.0` | `< 2.0` |
| **FS** | Fisher strand bias — are alt reads stuck on one strand? | `> 60` | `> 200` |
| **SOR** | another strand-bias measure (odds ratio) | `> 3.0` | `> 10.0` |
| **MQ** | root-mean-square mapping quality of the reads | `< 40` | — |
| **MQRankSum** | mapping-quality difference between ref and alt reads | `< −12.5` | — |
| **ReadPosRankSum** | are alt alleles biased to read ends (often artefacts)? | `< −8.0` | `< −20.0` |

Indels get **looser** strand/position thresholds because they are intrinsically
harder to align, which inflates these metrics even for real indels. Full reference:
[GATK annotations](https://gatk.broadinstitute.org/hc/en-us/articles/360035890471-Hard-filtering-germline-short-variants).

**Quick vocabulary:**
- **Read** – one short stretch of sequenced DNA. **Paired-end** = both ends read.
- **Reference / 3D7** – the standard genome we compare against.
- **Mapping / alignment** – placing a read at its genome position.
- **BAM** – compressed binary alignment file (`samtools` reads it).
- **Variant / SNP / indel** – a difference from the reference (single base / insertion-deletion).
- **VCF / GVCF** – the file formats that store variants / all-positions.
- **Codon** – 3 DNA bases coding one amino acid; resistance markers are named by codon (e.g. K76T = lysine→threonine at codon 76).
- **Coverage / depth** – reads overlapping a position.
- **Copy number** – how many copies of a gene the parasite has.

---

## Troubleshooting

| Symptom | Fix |
|---------|-----|
| `command not found: gatk` / `bwa` | Environment not active. A new terminal should show `(pfdr)`; else run the activate line from Part 1. |
| `ERROR: missing <file> - run step NN first` | You skipped a step or one failed. Run the named step; each script says what it needs. |
| Step 04 is slow | The ~3.5 GB human index is downloading (first time only). Let it finish. |
| `No space left on device` | Delete finished stages, e.g. `rm -rf results/01_fastqc_raw results/02_trim`, or use *Workshop mode* subsampling. |
| Report row says `coord_mismatch` | A guard, not a crash: the reference and annotation versions disagree at that codon. See `config.sh` to pin a PlasmoDB release. |
| An HTML report won't open | It's in the cloud VM — right-click ▸ *Open with Live Preview*, or download it and open locally. |

More in **[INSTRUCTORS.md](INSTRUCTORS.md)** (for facilitators).

---

## Running locally (optional, no Codespaces)

With [micromamba](https://mamba.readthedocs.io/en/latest/) or conda/mamba:

```bash
micromamba create -y -f environment.yml
micromamba activate pfdr
bash steps/00_reference.sh
# ...then the same commands as above
```

Needs internet (PlasmoDB + ENA + the hostile index), ~4 CPUs, ~16 GB RAM, a few GB
of disk.

---

## Repository layout

```
.devcontainer/      Codespaces setup (installs the pfdr environment) — don't delete
environment.yml     all tools + exact versions
config.sh           every setting, with comments
steps/              00–12, one script per stage (identical to the README commands)
scripts/            Python helpers (region builder, DR report, coverage plot)
resources/          dr_genes.tsv, dr_markers.tsv  (editable reference tables)
tests/              synthetic end-to-end self-test with known answers
run_all.sh          run every step in order
README.md           this course
INSTRUCTORS.md      facilitation guide: timing, class setup, prebuilds, troubleshooting
```

**Instructor self-test:** `bash tests/selftest.sh` builds a tiny synthetic genome
with planted mutations (*crt* K76T on the + strand, *dhfr* S108N on the − strand)
and a duplicated *mdr1*, runs the real step scripts, and checks the report and
copy-number estimate. Runs in under a minute, no downloads — a reliable smoke test.

---

## Going further

- **[MalariaGEN Pf7](https://www.malariagen.net/parasite/pf7)** — >20,000 *P.
  falciparum* genomes with curated resistance classifications; the reference
  resource for the field ([paper](https://doi.org/10.12688/wellcomeopenres.18681.1)).
- **[WHO Global Malaria Programme](https://www.who.int/teams/global-malaria-programme)**
  and the **[Malaria Threats Map](https://apps.who.int/malaria/maps/threats/)** —
  policy context and where markers are spreading.
- **[WWARN](https://www.wwarn.org/)** — antimalarial resistance surveillance network.
- **[Galaxy Training: variant analysis](https://training.galaxyproject.org/training-material/topics/variant-analysis/)**
  — a gentler, GUI-based take on the same concepts if you want more practice.
- **[GATK Best Practices](https://gatk.broadinstitute.org/hc/en-us/sections/360007226651-Best-Practices-Workflows)**
  — the canonical variant-calling workflows this course is based on.

> Reminder: this pipeline is for **training and research**, not clinical diagnosis.
