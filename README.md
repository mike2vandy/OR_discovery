# OR discovery workflow

## This workflow was used to discover, extract, and annotate Olfactor Receptor (OR) genes from turtle genomes in Cockrin et al 20xx. OR sequences are extracted from whole genome drafts and annotated to the level of intact, truncated and pseudogene sequences.
### The current workflow is written in snakemake 9.27 utilizing conda environments present in `env/`. The workflow is designed to perform best on a scheduled HPC with slurm, LSF, etc., but could run on a local desktop or server but may require resource re-configuration.

### The following prerequisites are required:
1. conda or mamba
2. snakemake v9
3. A local copy of DeepTMHMM.

### Installing snakemake: 
`conda create -n snakemake -c conda-forge snakemake=9`  

#### I ran this workflow with SLURM and snakemake-executor-plugins, which required.
```
conda activate snakemake
pip install snakemake-executor-plugin-slurm
```
#### However you may need to be adjust the install for your particular scheduler, i.e.:
`pip install snakemake-executor-plugin-<scheduler>`

### Installing DeepTMHMM:

1. Download a local academic version of DeepTMHMM from: `https://dtu.biolib.com/DeepTMHMM`. 

  - `unzip DeepTMHMM-Academic-License-v1.0.zip`

2. DeepTMHMM is CPU greedy. Once unzipped, you'll need to modify a few lines in `predict.py` to ensure pytorch wont reserve all cpus on an HPC node, but instead reserves the number of threads dictated by the snakemake rule deep_tm_hmm: 

```
#Update these two lines lines:

torch.set_num_threads(os.cpu_count())
torch.set_num_interop_threads(os.cpu_count() + 2)

#To 
n_threads = int(os.environ.get("DEEPTMHMM_THREADS", 1))

torch.set_num_threads(n_threads)
torch.set_num_interop_threads(n_threads + 2)
```

3. Update the `deep_tm_hmm_dir` variable in your config file to the full path of your local DeepTMHMM.
   -   `deep_tmhmm_dir : /path/to/DeepTMHMM-Academic-License-v1.0`

## Notes about a config file i.e. `turtle.yaml`:
### There are three customizable variables that can be modified for different uses:
```
input_table    : table.tst.csv 
deep_tmhmm_dir : /path/to//DeepTMHMM-Academic-License-v1.0
query_fasta    : query/intact.ORs.fas
```

  - input_table. this is a variable that holds the name of a csv file containing names of genome.fasta files you want to search, and short names that will be used as prefixes throughout the workflow. To query all turtle genomes in Cockrin et al. change to table.full.csv
  - deep_tmhmm_dir. See above
  - query_fasta. Set this variable to the path of a fasta file that will serve as the blast query. Currently intact.ORs.fas are known Sauropsid ORs. 

## Table notes:
### snakemake is controlled by the csv file assigned to `input_table`. This csv file should only contain two columns for the complete filename and a shortname.
```
genome,prefix
GCA_003846335.1_ASM384633v1_genomic.fna,cuomcc
GCA_003942145.1_ASM394214v1_genomic.fna,plameg
GCA_004028625.2_ASM402862v2_genomic.fna,cuoamb
```
Recognize the path is not needed and any genome can be included. Also recognize snakemake is expecting genome fasta files to be stored in the `genomes/` directory.

## Running the snakemake command:

The main command used was:
```
snakemake \
        -s OR_discovery.smk \
        --configfile turtles.yaml \ #any modified yaml file can be used. 
        --sdm conda \ 
        --executor slurm \
        --jobs 20 # number of jobs allowed to simultaneously run on a scheduler
```

 
## Replicating results in Cockrin et al 20xx:
### Once the prerequisites are installed:

1. Enter `genomes/` and run the `get_genomes.sh` bash script to download and unzip all queried turtle genomes.
   - `bash get_genomes.sh`
2. Return to main directory, create a job submission script for your scheduler and add the snakemake line.
   - Included is a `turtle_test.yaml` config file where `table.tst.csv` has a smaller subset of genomes for testing purposes. Run it like:
     ```
     snakemake \
        -s OR_discovery.smk \
        --configfile turtle_test.yaml 
        --sdm conda \ 
        --executor slurm \
        --jobs 20 \
     ```
   - You can include the `-n` flag to perform a dry run to ensure the workflow is performing properly before submitting.

   - To process all genomes:
     ```
     snakemake \
        -s OR_discovery.smk \
        --configfile turtles.yaml 
        --sdm conda \ 
        --executor slurm \
        --jobs 20 
     ``` 

## Output:
Results for each genome will be in `output/<prefix>/`. A completed sample should resembl:
```
output/tertri/
├── beds
│   ├── intact
│   │   └── tertri.intact_1.bed
│   ├── main
│   │   └── tertri.merged.bed
│   ├── pseudo
│   │   ├── tertri.pseudo_1.bed
│   │   └── tertri.pseudo_2.bed
│   └── truncated
│       ├── tertri.truncated_1.bed
│       └── tertri.truncated_2.bed
├── blast
│   └── tertri.blast.out
├── final
│   ├── tertri.complete_intact.fas
│   ├── tertri.complete_pseudo.bed
│   └── tertri.complete_truncated.bed
└── intact_fas
    ├── deep_tm_out
    │   ├── TMRs.gff3
    ├── tertri.intact_2.fas.transdecoder_dir
    │   ├── longest_orf.pep
    │   ├── longest_orf.cds
    │   ├── longest_orf.gff
    ├── tertri.complete_intact.fas
    ├── tertri.intact_1.fas
    ├── tertri.intact_2.fas
    ├── tertri.intact_2.fas.transdecoder.bed
    ├── tertri.intact_2.fas.transdecoder.cds
    ├── tertri.intact_2.fas.transdecoder.gff3
    └── tertri.intact_2.fas.transdecoder.pep
```
### The primary output is in `output/final` and should reflect:

  - <species_prefix>.complete_intact.fas: full length OR seqs with 7 transmembrane domains in fasta format and includes <prefix> in the fasta headers
  - <species_prefix>.copmplete_pseudo.bed: genome coordinates of OR genes that failed to meet intact parameters.
  - <species_prefix>.complete_truncated.bed: genome coordinates of OR genes that were too close to the end of a chromosome or assembly gap to assign intact or pseudogene status.

## Customization:



## General flowchart of the OR discory pipeline.

![](images/flowchart.png) 


To annotate ORs for you own purposes, store unique genomes in fasta format in `genomes`, modify `table.csv` to reflect genome names and desired prefixes, provide a query of OR aa sequences in fasta format to replace `query/intact.ORs.fas`
