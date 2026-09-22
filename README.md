# OR discovery workflow

## This workflow was used to discover, extract, and annotate in Cockrin et al 20xx. OR seqs are annotated to the level of intact, truncated and pseudogene sequences.
### workflow is written in snakemake 9.27 with conda environments present in `env/`

### Prerequisites

#### conda
#### snakemake 9
#### deeptmhmm
A local academic version of deeptmhmm is needed from `https://dtu.biolib.com/DeepTMHMM`. All depencies are in `envs/deeptmhmm.yaml. Unzip DeepTMHMM 
`unzip DeepTMHMM-Academic-License-v1.0.zip`

Update the path of DEEPTMHMM_DIR in config.  
The variable `DEEPTMHMM_DIR` in the snakemake script will need to be modified to reflect the full path of the deeptmhmm directory.  

update predict.py code with:

n_threads = int(os.environ.get("DEEPTMHMM_THREADS", 1))

torch.set_num_threads(n_threads)
torch.set_num_interop_threads(n_threads + 2)

pip executor

## Table notes

## config file
 
## Steps to run to replicate results in Cockrin et al 20xx.

1. Create a conda environment that contains snakemake 9.
   - `conda install mamba`
   - `mamba create -n snakemake -c bioconda -c conda-forge snakemake=9` 

1. Enter `genomes/` and run the `get_genomes.sh` bash script to download and unzip turtle genomes.
   - `bash get_genomes.sh`
2. Return to main directory, run the snakemake workflow
   - `snakemake -s OR_discovery.smk --cores 32 --use-conda`

## Results for each species will be in output/<species_prefix>/final/
  - <species_prefix>.complete_intact.fas: full length OR seqs with 7 transmembrane domains in fasta format and includes <species_prefix> in the header
  - <species_prefix>.copmplete_pseudo.bed: genome coordinates of OR genes that failed to meet intact parameters.
  - <species_prefix>.complete_truncated.bed: genome coordinates of OR genes that were too close to the end of a chromosome or assembly gap to assign intact or pseudogene status.

### General flowchart of the OR discory pipeline.

![](images/flowchart.png) 

To annotate ORs for you own purposes, store unique genomes in fasta format in `genomes`, modify `table.csv` to reflect genome names and desired prefixes, provide a query of OR aa sequences in fasta format to replace `query/intact.ORs.fas`
