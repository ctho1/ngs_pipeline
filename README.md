# NGS Pipeline
Customized bioinformatic pipeline for analyzing targeted NGS sequencing data and metagenomic sequencing experiments. Takes paired-end `.fastq.gz` files as input and generates various tabular output files and plots for visualization. The pipeline is adapted to the SLURM job scheduler for parallel processing of multiple samples. Requires 80 GB memory and adjustable number of CPUs.

# Arriba
Performs `STAR` alignment and downstream analysis with `Arriba` (https://github.com/suhrig/arriba) for gene fusion detection. Besides the human reference genome, reads are also aligned against ~12,000 RefSeq virus genomes (as well as ~4,500 human-infecting virus strains related to the RefSeq viruses) and viral integration sites into the human genome are also being analyzed.

# CNVkit
Both on- and off-target reads are used to generate tumor-only CNV profiles using `CNVkit` (https://github.com/etal/cnvkit). Output files include a genome-wide copy-number profile as well as detailed plots for each chromosome highlighting some tumor-relevant genes that are part of the NGS panel being currently used in our institute.

# Kraken
`Kraken2` (https://github.com/DerrickWood/kraken2) is being used for taxonomic classification of both off-target reads from targeted panel NGS data and direct metagenomics sequencing. In our hands, Kraken2 is a very sensitive, competetive and low-ressource tool for pathgen detection in diagnostic neuropathology. 
