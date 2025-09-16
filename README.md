# WDL-Somalier

This repository provides a WDL workflow for running [Somalier](https://github.com/brentp/somalier) to extract sample information from BAM files, relate samples, and check for identity, sex, and relationships using pedigree information.

## Features

- Batch extraction of Somalier files from BAMs
- Relate samples within and across batches
- Identity, sex, and relationship checks using pedigree data
- SLURM and Singularity support for scalable cluster execution
- Call caching to avoid redundant computation

## Input Files

- **Reference FASTA**: Path to GRCh38 reference genome
- **Sites VCF**: Path to variant sites for Somalier
- **Pedigree**: PED file for relationship checking
- **omeList**: List of sample "omes" to process
  - `ome_name`: is name of ome to organize by
  - `toExtractList`: TSV with `[ID, BAM, BAI]` for samples to extract
  - `extractedList`: TXT with paths to already extracted `.somalier` files

See [input.json](input.json) and [test_data/input.json](test_data/input.json) for examples.

## Usage

1. Prepare your input files as described above.
2. Run the workflow with miniwdl:

   ```sh
   miniwdl run somalier.wdl -v -d work -i input.json --cfg miniwdl.cfg

## Dependencies
# This repo assumes you have minwdl installed on HPC. If not, run the following commands
pip3 install miniwdl miniwdl-slurm
# add the path of miniwdl to your ~/.bashrc. 
which mindwdl
ls ~/.local/bin/miniwdl
