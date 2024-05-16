# WDL-Somalier

Example command  
miniwdl run somalier.wdl -v -d work -i input.json --cfg miniwdl.cfg

Input  
Add in GRCh38 reference for Somalier.reference  
omeList is a list of omes to process  
- ome_name is name of ome to organize by
- toExtractList is a tsv with columns [ID, path to BAM, path to index]. This file represents files that have not yet been extracted.
- extractedList is a txt with one .somalier file per row. This file represents already extracted samples.  
pedigree is a pedigree of the samples to use for relationship checking
