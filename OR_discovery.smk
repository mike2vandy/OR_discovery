import pandas as pd

SAMPLES = pd.read_csv('table.csv')
GENOMES = dict(zip(SAMPLES.prefix, SAMPLES.genome))

PREFIXES = SAMPLES.prefix.tolist()

rule all:
  input:
    expand('output/{prefix}/final/{prefix}.complete_truncated.bed', prefix = PREFIXES),
    expand('output/{prefix}/final/{prefix}.complete_pseudo.bed', prefix = PREFIXES),
    expand('output/{prefix}/final/{prefix}.complete_intact.fas', prefix = PREFIXES)
    
rule prepare_genomes:
  input:
    fas = lambda wildcards: f"genomes/{GENOMES[wildcards.prefix]}"
  output:
    db = "genome_dbs/{prefix}.ndb",
    fai = "genome_dbs/{prefix}.fai"
  params:
    db_type = "nucl",
    genome_prefix = "genome_dbs/{prefix}"
  conda: "env/or_env.yaml"
  threads: 1
  resources:
    mem_mb = 20000,
    time = 20
  shell:
    '''
       makeblastdb \
         -in {input.fas} \
         -out {params.genome_prefix} \
         -dbtype {params.db_type}

       samtools faidx \
         {input.fas} \
          -o {output.fai}
    '''

rule tblastn:
  input:
    db = "genome_dbs/{prefix}.ndb",
    query = "query/intact.ORs.fas"
  output:
    blast = "output/{prefix}/blast/{prefix}.blast.out",
  params:
    genome_prefix = "genome_dbs/{prefix}",
    evalue = '1e-10',
    fmt = '6'
  conda: "env/or_env.yaml"
  threads: 16
  resources:
    meem_mb = 40000,
    time = 1200
  shell:
    '''
      tblastn \
        -query {input.query} \
        -db {params.genome_prefix} \
        -evalue {params.evalue} \
        -outfmt {params.fmt} \
        -num_threads {threads} \
        -out {output.blast}
    '''

rule blast_to_bed:
  input:
    blast = "output/{prefix}/blast/{prefix}.blast.out"
  output:
    bed = "output/{prefix}/beds/main/{prefix}.merged.bed"
  conda: "env/or_env.yaml"
  threads: 1
  resources:
    mem_mb = 20000,
    time = 10
  shell:
    '''
      awk '{{if ($9 > $10) {{print $2"\t"$10"\t"$9}} \
           else {{print $2"\t"$9"\t"$10}}}}' {input.blast} \
        |bedtools sort \
        |bedtools merge > {output.bed}
    '''

rule sort_hits:
  input:
    fai = "genome_dbs/{prefix}.fai",
    bed = "output/{prefix}/beds/main/{prefix}.merged.bed"
  output:
    trunc = "output/{prefix}/beds/truncated/{prefix}.truncated_1.bed",
    pseud = "output/{prefix}/beds/pseudo/{prefix}.pseudo_1.bed", 
    intact = "output/{prefix}/beds/intact/{prefix}.intact_1.bed",
  threads: 1
  resources:
    mem_mb = 20000,
    time = 10
  run:    
    #define global variables
    genome_dict = {}
    flanking_seq = 150 
    truncated_file = open(output.trunc, 'w')
    pseudo_file = open(output.pseud, 'w')
    intact_file = open(output.intact, 'w') 
    
    #read in genome index file made by samtools faidx
    with open(input.fai, 'r') as f:
      for line in f:
        columns = line.strip().split()
        chrm, length = columns[0], columns[1]
        genome_dict[chrm] = int(length)
    
    #read through bed file, write out 3 files: truncated gene file, pseudogene file, prelim intact file
    
    with open(input.bed, 'r') as bed_file:
      for line in bed_file:
        columns = line.strip().split("\t")
        chrom = columns[0]
        start, end = int(columns[1]), int(columns[2])
        hit_length = end - start
        chr_len = genome_dict[chrom] #extract the chromosome length from the dictionary
        #conditional statement to ID truncated hits
        if start <= flanking_seq or end >= chr_len - flanking_seq:
          #write to truncated file
          print(f"{chrom}\t{start}\t{end}", file=truncated_file)
        #since these aren't truncated, we'll make decisions about them based on length:
        else:
          if hit_length <= 150:
            continue
          elif hit_length > 150 and hit_length <= 750:
            #write to pseudogene file
            print(f"{chrom}\t{start-flanking_seq}\t{end+flanking_seq}", file=pseudo_file)
          else:
            #write to an intact gene file
            print(f"{chrom}\t{start-flanking_seq}\t{end+flanking_seq}", file=intact_file)

    truncated_file.close()
    pseudo_file.close()
    intact_file.close()

rule extract_intact:
  input:
    genome = lambda wildcards: f"genomes/{GENOMES[wildcards.prefix]}",
    bed = "output/{prefix}/beds/intact/{prefix}.intact_1.bed"
  output:
    fasta = "output/{prefix}/intact_fas/{prefix}.intact_1.fas"
  conda: "env/or_env.yaml"
  threads: 1
  resources:
    mem_mb = 20000,
    time = 10
  shell:
    '''
      bedtools getfasta \
        -fi {input.genome} \
        -bed {input.bed} \
        -fo {output.fasta}
    '''
   
rule remove_truncated:
  input:
    fasta = "output/{prefix}/intact_fas/{prefix}.intact_1.fas"
  output:
    trunc = "output/{prefix}/beds/truncated/{prefix}.truncated_2.bed",
    intact = 'output/{prefix}/intact_fas/{prefix}.intact_2.fas'
  threads: 1
  resources:
    mem_mb = 20000,
    time = 10
  run:
    #files
    truncated_file = open(output.trunc, 'a')
    intact_file = open(output.intact, 'w')
  
    #load in intact fasta file, save seqs 
    seqs = {} 
    with open(input.fasta) as infile:
      for line in infile:
        line = line.strip()
        if line.startswith('>'):
          header = line
          seqs[header] = ''
        else:
          seqs[header] += line
    
    #check for assembly gaps at flanks, if present truncated 
    for header, seq in seqs.items():
      if seq.startswith('NNNNN') or seq.startswith('XXXXX') or seq.endswith('NNNNN') or seq.endswith('XXXXX'):
        coord = header.replace('>','').replace(':','\t').replace('-','\t')
        print(coord, file = truncated_file)
      else:
        print(f"{header}\n{seq}", file = intact_file)
    truncated_file.close()
    intact_file.close()

rule extract_orf:
  input:
    intact = 'output/{prefix}/intact_fas/{prefix}.intact_2.fas'
  output:
    aa = 'output/{prefix}/intact_fas/{prefix}.intact_2.fas.transdecoder_dir/longest_orfs.pep'
  params:
    output_dir = 'output/{prefix}/intact_fas/',
    aa_length = 250
  conda: 'env/or_env.yaml'
  threads: 1
  resources:
    mem_mb = 20000,
    time = 20
  shell:
    '''
      TransDecoder.LongOrfs \
        --complete_orfs_only \
        -m {params.aa_length} \
        -t {input.intact} \
        --output_dir {params.output_dir}
    ''' 

DEEPTMHMM_DIR = "/usr/local/usrapps/stern/mwvandew/DeepTMHMM-Academic-License-v1.0" 

rule deep_tm_hmm:
  input:
    pep = 'output/{prefix}/intact_fas/{prefix}.intact_2.fas.transdecoder_dir/longest_orfs.pep'
  output:
    tmm = directory("output/{prefix}/intact_fas/deep_tm_out")
  params:
    app = DEEPTMHMM_DIR 
  conda: 'env/deeptmhmm.yaml'
  threads: 32
  resources:
    gpu = 8,
    mem_mb = 60000,
    time = 1200
  shell:
    '''
      (
      cd {params.app}

      python predict.py \
        --fasta {workflow.basedir}/{input.pep} \
        --output-dir {workflow.basedir}/{output.tmm}
      )
    '''

rule parse_tm:
  input:
    tmm = rules.deep_tm_hmm.output.tmm,
    peps = 'output/{prefix}/intact_fas/{prefix}.intact_2.fas.transdecoder_dir/longest_orfs.pep', 
    prelim_intact = 'output/{prefix}/intact_fas/{prefix}.intact_2.fas' 
  output:
    intact = 'output/{prefix}/intact_fas/{prefix}.complete_intact.fas',
    pseudo = "output/{prefix}/beds/pseudo/{prefix}.pseudo_2.bed", 
  threads: 1
  resources:
    meme_mb = 20000,
    time = 20
  run:
    #load in transdecoder peps
    aas = {}
    with open(input.peps) as f:
      for line in f:
        line = line.strip()
        if line.startswith('>'):
          head = line.replace('>','').split()[0]
          aas[head] = ''
        else:
          aas[head] += line

    #load in all non-truncated intact seqs
    seqs = {}
    with open(input.prelim_intact) as f:
      for line in f:
        line = line.strip()
        if line.startswith('>'):
          head = line.replace('>','')
          seqs[head] = ''
        else:
          seqs[head] += line

    complete = {}
    pseudos = set()

    #parse gff from deep_tm_hmm
    with open(f"{input.tmm}/TMRs.gff3") as f:
      for line in f:
        line = line.strip()
        if "predicted" in line:
          fields = line.split()
          aa, tmds = fields[1], fields[-1]
          tmds = int(tmds)
          if tmds == 7: #save seqs with 7 TMDs
            complete[aa] = aas[aa]

    # if not in complete, write to pseudogene file as bed
    # going to account for peps < 250aa also
    with open(output.pseudo, 'w') as out:
      for i in seqs:
        key = i + '.p1'
        if key not in complete:
          seq = i.replace(':','\t').replace('-','\t')
          out.write(f"{seq}\n")

    #write complete orf seqs to intact
    with open(output.intact, 'w') as out2:
      for i, j in complete.items():
        out2.write(f">{i}\n{j}\n")
 
rule final_gather:
  input:
    trunc1 = "output/{prefix}/beds/truncated/{prefix}.truncated_1.bed",
    trunc2 = "output/{prefix}/beds/truncated/{prefix}.truncated_2.bed",
    pseud1 = "output/{prefix}/beds/pseudo/{prefix}.pseudo_1.bed",
    pseud2 = "output/{prefix}/beds/pseudo/{prefix}.pseudo_2.bed",
    intact = 'output/{prefix}/intact_fas/{prefix}.complete_intact.fas',
  output:
    final_trunc = 'output/{prefix}/final/{prefix}.complete_truncated.bed',
    final_pseud = 'output/{prefix}/final/{prefix}.complete_pseudo.bed',
    final_intact = 'output/{prefix}/final/{prefix}.complete_intact.fas',
  params:
    prefix = "{prefix}"
  threads: 1
  shell:
    '''
      cat {input.trunc1} {input.trunc2} > {output.final_trunc}
      cat {input.pseud1} {input.pseud2} > {output.final_pseud}

      sed 's/>/>{params.prefix}_/g' {input.intact} > {output.final_intact} 
    '''
