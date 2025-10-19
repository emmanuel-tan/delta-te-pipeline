#!/bin/bash

wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_30/gencode.v30.pc_transcripts.fa.gz
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_30/GRCh38.p12.genome.fa.gz
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_30/gencode.v30.primary_assembly.annotation.gtf.gz

gunzip gencode.v30.pc_transcripts.fa.gz
gunzip GRCh38.p12.genome.fa.gz
gunzip gencode.v30.primary_assembly.annotation.gtf.gz