# Instructor / facilitator guide

Everything you need to run this as a **2-hour hands-on course** for a class. Read
this once before the session; students only need `README.md`.

---

## 1. What the course teaches

By the end, a participant can take raw Illumina reads from a *P. falciparum*
sample and produce a **drug-resistance genotype**, understanding *why* each step
exists. Learning objectives are listed at the top of `README.md`.

The pipeline is deliberately **run by hand, step by step** (no Snakemake/Nextflow)
so students see and type every real command. The `steps/*.sh` scripts contain the
identical commands for reference and for the automated self-test.

---

## 2. Is 2 hours enough? (timing reality)

The **teaching content comfortably fills 2 hours** — the risk is the opposite:
**waiting on downloads/compute**. Two things dominate the clock:

| Wait | Size / time | Mitigation |
|------|-------------|------------|
| **hostile human index** (step 04) | ~3.5 GB download, once per Codespace | **Pre-cache it in a Codespaces *prebuild*** (§4) so it's already there, or have students start step 04 before a break |
| **Read download** (step 01) | the ENA FASTQs for ERR018972 | **Subsample** (§3), and/or prebuild-cache the reads |
| First Codespace build | ~3–5 min (building the `pfdr` env) | **Prebuilds** make this near-instant (§4) |

With prebuilds + subsampling, the hands-on time is dominated by *reading and
typing*, not waiting — which is what you want.

### Suggested 2-hour agenda

| Time | Segment | Steps |
|------|---------|-------|
| 0:00–0:15 | Welcome, the malaria drug-resistance problem, how Codespaces works, everyone launches their Codespace | — |
| 0:15–0:30 | Reference & reads; **FastQC** — read the data before touching it | 00, 01, 02 |
| 0:30–0:45 | Trimming (fastp) and host removal (why patient blood ≠ pure parasite) | 03, 04 |
| 0:45–1:05 | Mapping (`bwa mem -M`, read groups), duplicates, coverage | 05, 06, 07 |
| 1:05–1:30 | Variant calling (GVCF logic) + hard-filtering (and why not VQSR) | 08, 09 |
| 1:30–1:50 | The payoff: resistance report + mdr1 copy-number plot | 10, 11 |
| 1:50–2:00 | MultiQC, discussion, scaling to a cohort, Q&A | 12 |

If you're tight on time, **steps 04 and 08 are the two to start early / background**
(the index download and the caller are the slowest). Students can read the
explanation while those run.

---

## 3. Make it fast for a class: subsample

Full coverage is lovely but slow. For teaching, have everyone subsample at step 01:

```bash
SUBSAMPLE_PAIRS=2000000 bash steps/01_get_reads.sh
```

Everything downstream auto-detects the subsampled files. Trade-off: lower depth
weakens the **mdr1 copy-number** demo (step 11) and mixed-infection calls. Pick the
number from the sample's real depth — check it after step 07:

```bash
cat results/06_coverage/*.genome.mosdepth.summary.txt | awk '$1=="total"'
```

Aim to keep **≥ ~30× genome-wide** if you want the copy-number plot to be
convincing; drop lower only if you're short on time and will show the plot from
your own full run instead.

---

## 4. Set the class up on GitHub (recommended path)

**Goal:** every student gets their own independent, pre-warmed Codespace.

1. **Make this a template repo.** On the repo: *Settings ▸ General ▸ Template
   repository* ✔. Students click **Use this template ▸ Create a new repository**,
   then launch a Codespace from *their* copy. (Alternatively add them as
   collaborators, but template repos scale better and keep everyone isolated.)

2. **Turn on prebuilds** (this is the big time-saver). On the repo:
   *Settings ▸ Codespaces ▸ Set up prebuild ▸ Create*. GitHub builds the container
   (the whole `pfdr` environment) ahead of time, so students' Codespaces start in
   seconds instead of 3–5 minutes.
   Docs: <https://docs.github.com/en/codespaces/prebuilding-your-codespaces>

3. **Pre-cache the big downloads in the prebuild (optional but ideal).** To bake
   the ~3.5 GB human index (and optionally the reads) into the prebuilt image, add
   them to the create step. Append to `.devcontainer/setup.sh`:

   ```bash
   # Pre-fetch the hostile human index so step 04 is instant for students
   "$MM_BIN" run -n pfdr hostile index fetch --name human-t2t-hla --bowtie2 || true
   ```

   Rebuild the prebuild afterwards. (Caching the reads too: add a `wget` of the two
   FASTQs into `data/` here — but that enlarges the image; usually index-only is
   enough.)

4. **Costs / limits.** Free personal accounts get ~60 core-hours/month of
   Codespaces and 15 GB-month storage; this repo requests a 4-core machine, so the
   free tier is ~15 hours of active use. Tell students to **Stop** their Codespace
   when done (<https://github.com/codespaces>), and delete it after the course.
   GitHub Education / Teacher benefits raise these limits — worth applying for a
   class: <https://education.github.com/>.

5. **If `.devcontainer` didn't upload.** The web uploader silently drops hidden
   folders if files are hand-selected. Verify the repo shows a `.devcontainer/`
   folder with `devcontainer.json` and `setup.sh`. The safe way to upload is to
   **drag the whole `pf-dr-codespace` folder** onto the uploader in one go (that
   includes dotfiles), or use `git push`.

---

## 5. One-click launch badge

Add this to the top of `README.md` (replace `OWNER/REPO`):

```markdown
[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/OWNER/REPO)
```

`https://codespaces.new/OWNER/REPO` is the canonical one-click "create a Codespace"
URL. With a template repo, point students at **Use this template** instead so each
gets their own.

---

## 6. Common student problems (and fixes)

| Symptom | Cause | Fix |
|---------|-------|-----|
| `command not found: gatk` (etc.) | environment not active | new terminal should show `(pfdr)`; else `eval "$("$HOME/.local/bin/micromamba" shell hook -s bash)" && micromamba activate pfdr` |
| `ERROR: missing <file> - run step NN first` | steps run out of order, or a previous step failed | run the named step; each script states its prerequisite |
| Step 04 hangs/slow | 3.5 GB index downloading | expected on first run; prebuild-cache it (§4) |
| `No space left on device` | Codespace disk full (big index + BAMs) | `rm -rf results/0X_*` for stages you've finished, or subsample |
| Report shows `coord_mismatch` | GFF release ≠ FASTA release | it's a guard, not a crash; the self-check caught it. Pin both via `REF_FASTA_URL`/`REF_GFF_URL`/`PLASMODB_FALLBACK_RELEASE` in `config.sh` |
| GATK out-of-memory | small machine | lower `JAVA_MEM` in `config.sh`, or use a 4-core Codespace |
| Can't open an HTML report | it's in the cloud VM | right-click the file ▸ *Open with Live Preview*, or download it |

---

## 7. Check it works before the class

Run the **self-test** in a Codespace (or locally). It builds a synthetic genome
with planted mutations and runs the real step scripts — no downloads, < 1 minute:

```bash
bash tests/selftest.sh
```

Expect `Result: 10 passed, 0 failed`. This is your smoke test that the environment
and every script are healthy on the exact machine students will use.

---

## 8. Discussion prompts (fill the last 10 minutes)

- A patient's parasites show **mdr1 N86 (wild-type) and increased mdr1 copy
  number**. Which drug is most compromised, and why do single SNPs miss this?
- Why do we call with **ploidy 2** for a haploid parasite? What would a mixed
  infection look like in the VCF `AD` field?
- The sample is `crt K76T` positive. What does that predict for **chloroquine**,
  and why is chloroquine resistance now so widespread that K76T is almost a
  baseline in many regions?
- When would you switch from **hard-filtering** to **VQSR**? What do you need?
- Which steps change if the reads were **Oxford Nanopore** instead of Illumina?

---

## 9. Extension ideas (for a longer or follow-up session)

- Run a **second sample** and compare reports (change `SAMPLE` + the two URLs).
- Switch to **whole-genome** calling (`CALL_REGIONS=genome`) and discuss runtime.
- Add a gene to `resources/dr_markers.tsv` (e.g. a new *kelch13* codon) and re-run
  step 10 — shows how curated panels are maintained.
- Add **plasmepsin 2/3** copy-number plotting (clone step 11 with
  `--gene-id PF3D7_1408000`).
- Build a **joint-genotyped cohort** VCF with `GenomicsDBImport` (see README).

---

*This course is for training and research, not clinical diagnosis.*
