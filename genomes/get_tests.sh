#! /bin/bash
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/007/922/155/GCA_007922155.1_Mesoclemmys_tuberculata-1.0/GCA_007922155.1_Mesoclemmys_tuberculata-1.0_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/007/922/225/GCA_007922225.1_Emydura_subglobosa-1.0/GCA_007922225.1_Emydura_subglobosa-1.0_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/007/922/175/GCA_007922175.1_Pelusios_castaneus-1.0/GCA_007922175.1_Pelusios_castaneus-1.0_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/007/922/305/GCA_007922305.1_Dermatemys_mawii-1.0/GCA_007922305.1_Dermatemys_mawii-1.0_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/003/942/145/GCA_003942145.1_ASM394214v1/GCA_003942145.1_ASM394214v1_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/003/846/335/GCA_003846335.1_ASM384633v1/GCA_003846335.1_ASM384633v1_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/004/028/625/GCA_004028625.2_ASM402862v2/GCA_004028625.2_ASM402862v2_genomic.fna.gz
echo "finished downloading genomes"
echo "unzipping genome files..."
gunzip *.gz
echo "finished unzipping genomes"
