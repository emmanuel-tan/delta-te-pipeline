# A Pipeline for Performing Translational Efficiency Analysis From RNA-seq and Ribo-seq samples

> “End-to-end Ribo-seq + RNA-seq processing with contaminant removal, Salmon quant, tximeta/tximport gene counts, and deltaTE analysis.”

![main-diagram](assets/main-diagram.png)

# Introduction
To run translational efficiency analysis, various tools exist, but they are not often part of an integrated workflow specifically designed for this type of analysis. Pubicly available tools tend to exist standalone or may be part of a different pipeline. 

In order to address that, I built a pipeline that takes raw RNA-seq and Ribo-seq reads and processes them adequately to run deltaTE.

## Features at a glance
- **Dual-modality**: RIBO (SE) + RNA (PE)
- **Contaminant removal**: Bowtie2 against rRNA FASTA
- **Fast quantification**: Salmon (pseudoalignment)
- **Reusable indices**: **Skip index** build via `--salmon_index`
- **Reproducible**: implemented using Docker containers
- **Discoverable results**: gene counts, deltaTE gene tables, MultiQC reports

# Workflow 

![pipeline-diagram](assets/pipeline-diagram.png)

## Pipeline steps
1) Takes in single-end RPF and paired-end RNA
2) Removes contaminants (Bowtie2)
3) Quantifies with Salmon (optionally reusing a prebuilt)
4) Converts to Salmon transcript counts to gene-level counts (tximeta/tximport + custom tx2gene)
5) Computes differentially translated genes with deltaTE
6) Aggregates QC with MultiQC
7) Stores all results in `./results`

# Ways to run
## Testing and Demo Run

A demo run is completed automatically on this repo using Github Actions; the demo run can also be completed locally by pulling the repo, downloading the demo dataset, and running the pipeline using the profile `ci_demo`. 

*Nextflow and Docker should already be installed on your device*

**Download demo data and references from S3 bucket**
```bash
curl -fsSL -o "demo-v0.1.tar.zst" "https://delta-te-demo.s3.ap-southeast-1.amazonaws.com/public/demo-v0.1.tar.zst"
curl -fsSL -o "demo-reference-v0.1.tar.zst" "https://delta-te-demo.s3.ap-southeast-1.amazonaws.com/public/demo-reference-v0.1.tar.zst"
```

The downloaded files should be unzipped and placed into the `test/demo` directory.

**Pull the repo and change directory**

```bash
git clone https://github.com/emmanuel-tan/delta-te-pipeline.git
cd delta-te-pipeline
```

**Run the pipeline with the `ci_demo` profile**

```bash
nextflow run main.nf -profile ci_demo
# if running on a device with arm-based cpu architecture (such as an M1 Macbook), include the arm profile as well
nextflow run main.nf -profile ci_demo,arm
```

**Automated demo run**

As part of continuous integration set up, this repo will attempt a complete run of the pipeline using Github Actions when new pushes are made. It obtains demo data from an AWS S3 bucket containing three replicates of two samples for both RNA-seq and Ribo-seq data, with 1 million reads per sample. It runs quantification on a prebuilt Salmon index on a randomly sampled one-third of the human reference transcriptome due to size constraints.

After completing the pipeline, the Github Actions workflow also publishes the reports on Github Pages which can be viewed [here](https://emmanuel-tan.github.io/delta-te-pipeline/).

This serves as part of a testing framework but also as a demo run of the pipeline. It takes about 20 minutes to complete this run on the free Github Actions Ubuntu runner, which has the following specifications [as stated by Github](https://docs.github.com/en/actions/reference/runners/github-hosted-runners#standard-github-hosted-runners-for-public-repositories). 

| Virtual Machine | Processor (CPU) | Memory (RAM) | Storage (SSD) | Architecture | Workflow label |
| ----- | ----- | ----- | ----- | ----- | ----- |
| Linux | 4 | 16 GB | 14 GB | x64 | ubuntu-latest, ubuntu-24.04, ubuntu-22.04 |

## Standard run 

To run this on your own, the following parameters are expected. 

| Parameter | Expected input |
| --- | --- |
| `--adapters`              | path to `fasta` file containing adapter sequences |
| `--abundantReference`     | path to `fasta` file containing sequences of contaminants to filter |
| `--riboseq_samplesheet`   | path to `csv` file containing Ribo-seq data samplesheet ([example](./test/demo/RIBO-samplesheet.csv)) |
| `--rnaseq_samplesheet`    | path to `csv` file containing RNA-seq data samplesheet ([example](./test/demo/RNA-samplesheet.csv)) |
| `--deltate_samplesheet`   | path to `tsv` file containing deltaTE samplesheet ([example](./test/demo/deltaTE-Samplesheet.tsv)) |
| `--humanReferenceTx`      | path to `fasta` reference transcriptome |
| `--humanReferenceGenome`  | path to `fasta` reference genome |
| `--humanReferenceAnnot`   | path to `gtf` reference annotations |
| `--salmon_index`          | path to pre-built Salmon index (optional) |

Note: While this pipeline can build a Salmon index, it is suggested to run the pipeline with a **pre-built Salmon index** on devices with low RAM available and to save time. 

## Outputs
All output of each step is saved to the `./results` directory. 

This includes
- `results/salmon/` — quant.sf, lib_format_counts.json, salmon_quant.log
- `results/tximport/` — gene-level count matrices (RIBO/RNA)
- `results/deltaTE/` — gene tables/plots
- `results/qc/` - per-sample QC and summarized MulitQC report (.html, .zip)


# Troubleshooting / FAQ

# Tools used

<details><summary>Nextflow</summary>
Di Tommaso, P., Chatzou, M., Floden, E. W., Barja, P. P., Palumbo, E., & Notredame, C. (2017). Nextflow enables reproducible computational workflows. <i>Nature biotechnology, 35</i>(4), 316-319.

- https://github.com/nextflow-io/nextflow
</details>

<details><summary>deltaTE</summary>
Chothani, S., Adami, E., Ouyang, J. F., Viswanathan, S., Hubner, N., Cook, S. A., ... & Rackham, O. J. (2019). deltaTE: Detection of translationally regulated genes by integrative analysis of Ribo‐seq and RNA‐seq data. <i>Current protocols in molecular biology, 129</i>(1), e108.

- https://github.com/SGDDNB/translational_regulation
</details>

<details><summary>Trimmomatic</summary>
Bolger, A. M., Lohse, M., & Usadel, B. (2014). Trimmomatic: a flexible trimmer for Illumina sequence data. <i>Bioinformatics, 30</i>(15), 2114-2120.

- https://github.com/usadellab/Trimmomatic
</details>

<details><summary>Bowtie2</summary>
Langmead B, Salzberg S. Fast gapped-read alignment with Bowtie 2. Nature Methods. 2012, 9:357-359.

- https://github.com/BenLangmead/bowtie2
</details>

<details><summary>Salmon</summary>
Patro, R., Duggal, G., Love, M. I., Irizarry, R. A., & Kingsford, C. (2017). Salmon provides fast and bias-aware quantification of transcript expression. <i>Nature methods, 14</i>(4), 417-419.

- https://combine-lab.github.io/salmon/
</details>

<details><summary>tximport</summary>
- https://github.com/thelovelab/tximport
</details>

<details><summary>FastQC</summary>
- https://github.com/s-andrews/FastQC
</details>

<details><summary>MultiQC</summary>
Ewels, P., Magnusson, M., Lundin, S., & Käller, M. (2016). MultiQC: summarize analysis results for multiple tools and samples in a single report. <i>Bioinformatics, 32</i>(19), 3047-3048.

- https://github.com/MultiQC/MultiQC
</details>


# Roadmap

While this project is functional and complete in its current form, it is not yet deployment-ready. The following items outline future improvements and enhancements planned to bring it closer to production, research, or commercial standards:

- [ ] Improved parameter management and finetuning (ie. CPU and RAM usage, custom directory to save results)
- [ ] Pre-checks to catch early formatting fails and pre-pulls of Docker images
- [ ] Profiles for running on HPCs or cloud executors

# Contact
If you have any feedback, feel free to contact me at emmanueltan2000@gmail.com, or send me a DM on [Linkedin](https://www.linkedin.com/in/emmanuel-tan-0b89051b3/).