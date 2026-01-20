# WDL-Somalier

Example command  
miniwdl run somalier.wdl -v -d work -i input.json --cfg miniwdl.cfg

Input  
Add in GRCh38 reference for Somalier.reference  
pedigree is a pedigree of the samples to use for relationship checking  
omeList is a list of omes to process  
- ome_name is name of ome to organize by
- toExtractList is a tsv with columns [ID, path to BAM, path to index]. This file represents files that have not yet been extracted.
- extractedList is a txt with one .somalier file per row. This file represents already extracted samples.  

## Features

### Identity Checking
The workflow performs two types of identity checks:
1. **Within-ome checking** (CheckIdentical): Verifies sample identity within each sequencing type (e.g., PacBio, RNA-seq, SR-DNA)
2. **Cross-ome checking** (CheckIdenticalAcrossOmes): Verifies sample identity across all sequencing types to ensure the same sample ID maintains consistent genetic identity regardless of sequencing platform

Both checks flag samples where relatedness falls below the identity threshold (default: 0.95).

