#!/usr/bin/env bash
# =====================================================================
#  run_all.sh
#  One script that runs the whole workshop exactly as written in YOUR
#  edited README (the attachment), steps 00 -> 12, end to end.
#
#  Run it from the repository root (same folder as the steps/ and
#  scripts/ directories), inside the Codespace with the pfdr env active:
#
#      bash run_my_pipeline.sh
#
#  Notes on the few tiny robustness tweaks vs your pasted commands
#  (each is marked [tweak] below so you can see exactly what changed):
#    - mkdir -p instead of mkdir            (so re-runs don't error)
#    - rm -f the .dict before CreateSequenceDictionary (it refuses to
#      overwrite an existing dictionary)
#    - gunzip -f Pf3D7.gff.gz   (unambiguous form of your `gunzip Pf3D7.gff`)
#  Nothing else is changed: same URLs, same flags, same paths.
# =====================================================================
set -euo pipefail

# run from the directory this script lives in (should be the repo root)
cd "$(dirname "$0")"

# ---- settings --------------------------------------------------------
S=ERR018972            # sample name
T=4                    # threads
SUBSAMPLE_PAIRS=2000000   # your step 03 reads data/${S}.sub_* , so we subsample
export SUBSAMPLE_PAIRS

banner() { echo; echo "############################################################"; echo "#  $*"; echo "############################################################"; }

# =====================================================================
banner "STEP 00  Reference genome + indexes + DR regions"
# =====================================================================
mkdir -p ref            # [tweak] -p
cd ref
wget -O Pf3D7.fasta "https://zenodo.org/records/18681980/files/PlasmoDB-54-Pfalciparum3D7-Genome.fasta?download=1"
wget -O Pf3D7.gff.gz  "https://zenodo.org/records/18681980/files/PlasmoDB-55_Pfalciparum3D7.gff.gz?download=1"
gunzip -f Pf3D7.gff.gz   # [tweak] -> Pf3D7.gff  (your line was `gunzip Pf3D7.gff`)

# Build indexes
samtools faidx Pf3D7.fasta
bwa index Pf3D7.fasta
rm -f Pf3D7.dict         # [tweak] CreateSequenceDictionary won't overwrite
gatk CreateSequenceDictionary -R Pf3D7.fasta

# Get DR BED + CDS model (run from repo root)
cd ../
python ./scripts/build_dr_regions.py \
      --fasta ./ref/Pf3D7.fasta \
      --gff ./ref/Pf3D7.gff \
      --genes ./resources/dr_genes.tsv \
      --out-bed ./resources/dr_genes.bed \
      --out-cds ./resources/dr_genes.cds.json

# =====================================================================
banner "STEP 01  Download reads (subsampled to ${SUBSAMPLE_PAIRS} pairs)"
# =====================================================================
# downloads the full FASTQs AND writes data/${S}.sub_* that step 03 uses
SUBSAMPLE_PAIRS=${SUBSAMPLE_PAIRS} bash steps/01_get_reads.sh

# =====================================================================
banner "STEP 02  FastQC on raw reads"
# =====================================================================
mkdir -p results/01_fastqc_raw
fastqc -t $T --outdir results/01_fastqc_raw \
    data/${S}_1.fastq.gz data/${S}_2.fastq.gz

# =====================================================================
banner "STEP 03  Trim with fastp, then re-FastQC"
# =====================================================================
mkdir -p results/02_trim
fastp \
    --in1 data/${S}.sub_1.fastq.gz --in2 data/${S}.sub_2.fastq.gz \
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

# =====================================================================
banner "STEP 04  Remove human reads with hostile"
# =====================================================================
mkdir -p results/03_dehost
# hostile writes outputs into the CURRENT dir (no --out-dir), then we move them.
# First run downloads the ~3.5 GB human-t2t-hla index and caches it.
hostile clean \
    --fastq1 results/02_trim/${S}.trim_1.fastq.gz \
    --fastq2 results/02_trim/${S}.trim_2.fastq.gz \
    --aligner bowtie2 \
    --index human-t2t-hla \
    --threads $T

mv ./${S}.trim_1.clean_1.fastq.gz ./results/03_dehost/${S}.clean_1.fastq.gz
mv ./${S}.trim_2.clean_2.fastq.gz ./results/03_dehost/${S}.clean_2.fastq.gz

# =====================================================================
banner "STEP 05  Map with bwa mem -M  ->  sorted BAM"
# =====================================================================
mkdir -p results/04_map
bwa mem -M -t $T \
    -R "@RG\tID:${S}\tSM:${S}\tLB:${S}_lib1\tPL:ILLUMINA" \
    ref/Pf3D7.fasta \
    results/03_dehost/${S}.clean_1.fastq.gz \
    results/03_dehost/${S}.clean_2.fastq.gz \
  | samtools sort -@ $T -o results/04_map/${S}.sorted.bam -
samtools index results/04_map/${S}.sorted.bam
samtools flagstat results/04_map/${S}.sorted.bam

# =====================================================================
banner "STEP 06  MarkDuplicates"
# =====================================================================
mkdir -p results/05_markdup
gatk MarkDuplicates \
    -I results/04_map/${S}.sorted.bam \
    -O results/05_markdup/${S}.markdup.bam \
    -M results/05_markdup/${S}.dup_metrics.txt \
    --CREATE_INDEX true

# =====================================================================
banner "STEP 07  Coverage with mosdepth"
# =====================================================================
mkdir -p results/06_coverage
mosdepth --no-per-base --fast-mode -t $T --by 500 \
    results/06_coverage/${S}.genome \
    results/05_markdup/${S}.markdup.bam
mosdepth --no-per-base -t $T --by ./resources/dr_genes.bed \
    results/06_coverage/${S}.dr \
    results/05_markdup/${S}.markdup.bam
echo "--- per-DR-gene mean depth ---"
zcat results/06_coverage/${S}.dr.regions.bed.gz
echo "--- genome-wide average ---"
awk '$1=="total"' results/06_coverage/${S}.genome.mosdepth.summary.txt

# =====================================================================
banner "STEP 08  Call variants (HaplotypeCaller GVCF -> GenotypeGVCFs)"
# =====================================================================
mkdir -p results/07_variants
gatk HaplotypeCaller \
    -R ref/Pf3D7.fasta \
    -I results/05_markdup/${S}.markdup.bam \
    -O results/07_variants/${S}.g.vcf.gz \
    -ERC GVCF \
    --sample-ploidy 2 \
    -L ./resources/dr_genes.bed --interval-padding 1000
gatk GenotypeGVCFs \
    -R ref/Pf3D7.fasta \
    -V results/07_variants/${S}.g.vcf.gz \
    -O results/07_variants/${S}.vcf.gz
echo "raw variant count:"
zcat results/07_variants/${S}.vcf.gz | grep -vc '^#' || true

# =====================================================================
banner "STEP 09  Hard-filter (split SNP/indel, filter, merge)"
# =====================================================================
mkdir -p results/08_filter
gatk SelectVariants -R ref/Pf3D7.fasta -V results/07_variants/${S}.vcf.gz \
    --select-type-to-include SNP   -O results/08_filter/${S}.snps.raw.vcf.gz
gatk SelectVariants -R ref/Pf3D7.fasta -V results/07_variants/${S}.vcf.gz \
    --select-type-to-include INDEL -O results/08_filter/${S}.indels.raw.vcf.gz

gatk VariantFiltration -R ref/Pf3D7.fasta \
    -V results/08_filter/${S}.snps.raw.vcf.gz \
    -O results/08_filter/${S}.snps.filt.vcf.gz \
    --filter-expression "QD < 2.0"              --filter-name "QD2" \
    --filter-expression "FS > 60.0"             --filter-name "FS60" \
    --filter-expression "MQ < 40.0"             --filter-name "MQ40" \
    --filter-expression "MQRankSum < -12.5"     --filter-name "MQRankSum" \
    --filter-expression "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum" \
    --filter-expression "SOR > 3.0"             --filter-name "SOR3"

gatk VariantFiltration -R ref/Pf3D7.fasta \
    -V results/08_filter/${S}.indels.raw.vcf.gz \
    -O results/08_filter/${S}.indels.filt.vcf.gz \
    --filter-expression "QD < 2.0"               --filter-name "QD2" \
    --filter-expression "FS > 200.0"             --filter-name "FS200" \
    --filter-expression "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum" \
    --filter-expression "SOR > 10.0"             --filter-name "SOR10"

gatk MergeVcfs \
    -I results/08_filter/${S}.snps.filt.vcf.gz \
    -I results/08_filter/${S}.indels.filt.vcf.gz \
    -O results/08_filter/${S}.filtered.vcf.gz
echo "--- FILTER breakdown ---"
zcat results/08_filter/${S}.filtered.vcf.gz | grep -v '^#' | cut -f7 | sort | uniq -c || true

# =====================================================================
banner "STEP 10  Drug-resistance report"
# =====================================================================
mkdir -p results/09_report
python scripts/report_dr.py \
    --fasta ref/Pf3D7.fasta \
    --cds ./resources/dr_genes.cds.json \
    --markers resources/dr_markers.tsv \
    --vcf results/08_filter/${S}.filtered.vcf.gz \
    --coverage results/06_coverage/${S}.dr.regions.bed.gz \
    --sample $S \
    --outdir results/09_report \
    --min-dp 5 \
    --mixed-min-af 0.10

# =====================================================================
banner "STEP 11  mdr1 coverage + copy-number plot"
# =====================================================================
mkdir -p results/10_plots
python scripts/plot_coverage.py \
    --windows results/06_coverage/${S}.genome.regions.bed.gz \
    --summary results/06_coverage/${S}.genome.mosdepth.summary.txt \
    --dr-regions results/06_coverage/${S}.dr.regions.bed.gz \
    --cds ./resources/dr_genes.cds.json \
    --sample $S \
    --gene mdr1 --gene-id PF3D7_0523000 \
    --out results/10_plots/${S}.mdr1_coverage.png

# =====================================================================
banner "STEP 12  MultiQC"
# =====================================================================
mkdir -p results/11_multiqc
multiqc -f -o results/11_multiqc -n ${S}_multiqc.html results

echo
echo "DONE. Key outputs:"
echo "  results/09_report/${S}.dr_variants.tsv"
echo "  results/10_plots/${S}.mdr1_coverage.png"
echo "  results/11_multiqc/${S}_multiqc.html"
