# ERICLoD

**ERICLoD** is a synthetic NGS data generator for validating the **Limit of Detection (LoD)** of clinically relevant somatic variants, aligned with the **2024 ERIC recommendations** for TP53 mutation analysis in CLL. This release includes full support for **TP53**. Support for other ERIC-relevant genes is under active development.

ERICLoD simulates sequencing data at user-defined VAF levels, read depths, and library prep parameters, allowing labs to benchmark their variant calling pipelines and determine their true analytical sensitivity under ERIC-aligned standards.

---

## 🚀 Key Features

- Simulates TP53 mutations at arbitrary VAFs (e.g., 1%, 5%, 10%)
- Generates realistic FASTQ/SAM output for WT and mutant alleles
- Supports paired-end, insert size, platform-specific ART simulation
- Built-in truth table and simulation summary output
- `--check-ref-match` mode for verifying FASTA-VCF compatibility
- Reproducible simulations with seed control
- Fully modular and CLI-driven, ready for integration into bioinformatics workflows

---

## 🧪 Why ERICLoD? The Rationale

The 2024 ERIC recommendations have redefined best practices in TP53 mutation analysis for CLL:

- Sanger sequencing is no longer sufficient
- All TP53 mutations must be reported, **regardless of VAF**
- Labs must determine and document their **LoD** through internal validation
- Validation must account for read length, coverage, and prep method
- Alignment accuracy to the correct reference genome must be verified

**ERICLoD enables labs to meet these demands** by simulating controlled synthetic datasets that reflect real-world sequencing challenges — with known ground truth.

---

## 📊 ERICLoD Capabilities vs. ERIC Recommendations

| **ERIC Recommendation**                                                                 | **ERICLoD Feature**                                                                                  |
|------------------------------------------------------------------------------------------|------------------------------------------------------------------------------------------------------|
| LoD must be validated with known inputs                                                 | Simulates known TP53 mutations at VAFs from 0.01 to 1.0                                            |
| Report all mutations regardless of VAF if read support is sufficient                    | Adjustable coverage and VAF, WT and mutant read simulation                                         |
| Bioinformatics workflows must be validated internally                                   | Generates FASTQs + truth sets for pipeline benchmarking                                            |
| Reference genome alignment accuracy must be confirmed                                   | `--check-ref-match` flags reference/ref-allele mismatches                                           |
| Validation must consider prep method, read length, and fragment size                    | Full control over ART settings: paired-end, insert size, platform                                  |
| Documentation required for audit/accreditation (e.g., ISO 15189)                        | Outputs: truth_table.tsv, simulation_summary.tsv, mismatch logs, zipped final output               |

---

## 🔧 Customizability

Although the current default dataset is based on **ClinVar**, labs can:

- Provide their own internal/private VCF files
- Simulate lab-specific variants
- Reuse internal variant panels for LoD assessment

ClinVar is used here **only because it is publicly accessible**. The engine works with any valid VCF input.

---

## 🖥️ Installation (Conda)

Create a conda environment and install dependencies with:

```bash
conda env create -f ericlod_env.yaml
conda activate ericlod
```

### `ericlod_env.yaml`:
```yaml
name: ericlod
channels:
  - conda-forge
  - bioconda
  - defaults
dependencies:
  - python=3.10
  - pyfaidx
  - art=2016.06.05
  - samtools
  - gzip
```

---

## 📂 Output Overview

ERICLoD outputs a fully packaged directory:

```
my_simulation/
├── region_1/
│   ├── original.fa
│   ├── mutated.fa
│   ├── wt1.fq, wt2.fq
│   ├── mut1.fq, mut2.fq
│   └── *.sam, *.aln
├── truth_table.tsv
├── simulation_summary.tsv
├── ref_mismatch_summary.tsv (if --check-ref-match enabled)
└── my_simulation.zip
```

---

## ⚖️ Ethical Use

### Disclaimer:
> **ERICLoD is not affiliated with or endorsed by the European Research Initiative on CLL (ERIC).** The tool is designed to align with the public recommendations published by ERIC for TP53 mutation analysis and LoD validation.

Use of this tool should always follow institutional, clinical, and regulatory guidelines. It is intended for research, validation, and quality assurance purposes.

---

## 📌 Citation

If you use ERICLoD in your research or validation workflow, please cite the following guideline document that forms the basis for this tool:

> Rossi D, Zenz T, Malcikova J, et al. **Revised ERIC recommendations for TP53 mutation analysis in chronic lymphocytic leukemia**. *Leukemia*. 2024. https://doi.org/10.1038/s41375-024-02135-3

These recommendations serve as the scientific foundation and clinical rationale for the existence and design of ERICLoD.

---

## 🧭 Example Usage

```bash
python3 ericlod.py \
  --vcf clinvar_tp53.vcf \
  --ref GRCh38.fa \
  --vaf 0.05 \
  --paired \
  --rules missense_variant=5,frameshift_variant=2 \
  --check-ref-match \
  --coverage 1000 \
  --output my_simulation
```

---

## 📌 Future Work

- Multi-gene support (ATM, NOTCH1, SF3B1)
- Graph-based mutation simulation
- LoD curve generation
- BAM/CRAM post-processing module

---

**ERICLoD** helps turn ERIC's recommendations into real-world diagnostics. You define the variant burden. We simulate the truth. You validate your precision.

---

🧬 *Built for translational research, clinical NGS, and high-stakes validation.*


