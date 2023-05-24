## Genome-Guided Transcriptome Assembly
The first step is to align the male reads to the genome using tophat after indexing the genome using Bowtie2:
```
module load tophat
module load bowtie2/2.4.4

srun bowtie2-build Artemia_sinica_genome_29_12_2021.fasta genome_w

tophat -p 40 \
    -o tophat_all \
    genome_w \
 398_1.fastq \
 398_2.fastq

```
The bam file from tophat is sorted using samtools and then used as input to Trinity:
```
module load java

module load samtools/1.10

module load jellyfish/2.3.0

module load bowtie2/2.4.4

module load python/3.8.5

module load salmon/0.13.1

samtools sort /nfs/scistore03/vicosgrp/melkrewi/genome_assembly_december_2021/774.genome_guided_transcriptome_assembly/tophat_all/accepted_hits.bam -o rnaseq.coordSorted.bam

/nfs/scistore03/vicosgrp/melkrewi/Schistosome_Project/10.transcriptome_assembly/data/Trinity/trinityrnaseq-v2.11.0/Trinity --genome_guided_bam rnaseq.coordSorted.bam --genome_guided_max_intron 10000 --max_memory 99G --CPU 20
```
## De-novo Transcriptome Assembly:
```
module load java

module load samtools/1.10

module load jellyfish/2.3.0

module load bowtie2/2.4.4

module load python/3.8.5

module load salmon/0.13.1

#run commands on SLURM's srun

/nfs/scistore03/vicosgrp/melkrewi/Schistosome_Project/10.transcriptome_assembly/data/Trinity/trinityrnaseq-v2.11.0/Trinity --seqType fq --left 39869_TTAGGC_C8W3EANXX_3_20160523B_20160523.1.fastq,39870_AGTTCC_C8W3EANXX_3_20160523B_20160523.1.fastq --right 39869_TTAGGC_C8W3EANXX_3_20160523B_20160523.2.fastq,39870_AGTTCC_C8W3EANXX_3_20160523B_20160523.2.fastq --CPU 20 --max_memory 99G
```

## Evigene:

```
module load java

module load blast

module load exonerate

module load cdhit

cat Trinity.fasta Trinity-GG.fasta > denovo_and_guided.fasta
perl /nfs/scistore03/vicosgrp/melkrewi/genome_assembly_december_2021/774.genome_guided_transcriptome_assembly/evigene/scripts/prot/tr2aacds.pl -cdnaseq denovo_and_guided.fasta
```

## Expression Analysis:

The expression analysis was perferformed using the cds from evigene, (only the isoform t1 was kept and with cds > 500 bp):
```
../faFilter -minSize=500 Artemia_sinica_transcriptome.cds Artemia_sinica_transcriptome_500bp.cds
```
## Kallisto:
```
module load kallisto

srun kallisto index -i transcripts_sinica.idx Artemia_sinica_transcriptome_500bp.cds

srun kallisto quant -t 40 -i transcripts_sinica.idx -o f1 -b 100 ./females/39871_ATTCCT_C8W3EANXX_4_20160523B_20160523.1.fastq ./females/39871_ATTCCT_C8W3EANXX_4_20160523B_20160523.2.fastq

srun kallisto quant -t 40 -i transcripts_sinica.idx -o f2 -b 100 ./females/39872_ATGTCA_C8W3EANXX_4_20160523B_20160523.1.fastq ./females/39872_ATGTCA_C8W3EANXX_4_20160523B_20160523.2.fastq

srun kallisto quant -t 40 -i transcripts_sinica.idx -o f3 -b 100 ./females/A.sinica.F.head.RNA.69_1.fq ./females/A.sinica.F.head.RNA.69_2.fq

srun kallisto quant -t 40 -i transcripts_sinica.idx -o f4 -b 100 ./females/A.sinica.F.head.RNA.70_1.fq ./females/A.sinica.F.head.RNA.70_2.fq

srun kallisto quant -t 40 -i transcripts_sinica.idx -o f5 -b 100 ./females/A.sinica.F.head.RNA.97_1.fq ./females/A.sinica.F.head.RNA.97_2.fq

srun kallisto quant -t 40 -i transcripts_sinica.idx -o f6 -b 100 ./females/A.sinica.F.head.RNA.98_1.fq ./females/A.sinica.F.head.RNA.98_2.fq

srun kallisto quant -t 40 -i transcripts_sinica.idx -o f7 -b 100 ./females/A.sinica.F.ovaries.RNA.40771_1.fq ./females/A.sinica.F.ovaries.RNA.40771_2.fq

srun kallisto quant -t 40 -i transcripts_sinica.idx -o f8 -b 100 ./females/A.sinica.F.ovaries.RNA.40772_1.fq ./females/A.sinica.F.ovaries.RNA.40772_2.fq

srun kallisto quant -t 40 -i transcripts_sinica.idx -o f9 -b 100 ./females/A.sinica.F.ovaries.RNA.90_1.fq ./females/A.sinica.F.ovaries.RNA.90_2.fq

srun kallisto quant -t 40 -i transcripts_sinica.idx -o f10 -b 100 ./females/A.sinica.F.ovaries.RNA.99_1.fq ./females/A.sinica.F.ovaries.RNA.99_2.fq

srun kallisto quant -t 40 -i transcripts_sinica.idx -o f11 -b 100 ./females/A.sinica.F.thorax.RNA.45052_1.fq ./females/A.sinica.F.thorax.RNA.45052_2.fq

srun kallisto quant -t 40 -i transcripts_sinica.idx -o f12 -b 100 ./females/A.sinica.F.thorax.RNA.45053_1.fq ./females/A.sinica.F.thorax.RNA.45053_2.fq

srun kallisto quant -t 40 -i transcripts_sinica.idx -o m1 -b 100 ./males/39869_TTAGGC_C8W3EANXX_3_20160523B_20160523.1.fastq ./males/39869_TTAGGC_C8W3EANXX_3_20160523B_20160523.2.fastq

srun kallisto quant -t 40 -i transcripts_sinica.idx -o m2 -b 100 ./males/39870_AGTTCC_C8W3EANXX_3_20160523B_20160523.1.fastq ./males/39870_AGTTCC_C8W3EANXX_3_20160523B_20160523.2.fastq

srun kallisto quant -t 40 -i transcripts_sinica.idx -o m3 -b 100 ./males/39877_ACAGTG_C8W3EANXX_4_20160523B_20160523.1.fastq ./males/39877_ACAGTG_C8W3EANXX_4_20160523B_20160523.2.fastq

srun kallisto quant -t 40 -i transcripts_sinica.idx -o m4 -b 100 ./males/39878_GGCTAC_C8W3EANXX_4_20160523B_20160523.1.fastq ./males/39878_GGCTAC_C8W3EANXX_4_20160523B_20160523.2.fastq

srun kallisto quant -t 40 -i transcripts_sinica.idx -o m5 -b 100 ./males/39879_TTAGGC_C8W3EANXX_4_20160523B_20160523.1.fastq ./males/39879_TTAGGC_C8W3EANXX_4_20160523B_20160523.2.fastq

srun kallisto quant -t 40 -i transcripts_sinica.idx -o m6 -b 100 ./males/39880_1.fq ./males/39880_2.fq

srun kallisto quant -t 40 -i transcripts_sinica.idx -o m7 -b 100 ./males/39895_1.fq ./males/39895_2.fq

srun kallisto quant -t 40 -i transcripts_sinica.idx -o m8 -b 100 ./males/39896_1.fq ./males/39896_2.fq

srun kallisto quant -t 40 -i transcripts_sinica.idx -o m9 -b 100 ./males/39901_1.fq ./males/39901_2.fq

srun kallisto quant -t 40 -i transcripts_sinica.idx -o m10 -b 100 ./males/39902_1.fq ./males/39902_2.fq

srun kallisto quant -t 40 -i transcripts_sinica.idx -o m11 -b 100 ./males/40767_1.fq ./males/40767_2.fq

srun kallisto quant -t 40 -i transcripts_sinica.idx -o m12 -b 100 ./males/40768_1.fq ./males/40768_2.fq
```
## normalization:
```
library(dplyr)
exp<-read.table("Expression_sinica.txt", head=T, sep=",")

expf<-exp[,-1]

rownames(expf)<-exp[,1]
par(mfrow=c(2,1))
par(mar=c(3,3,0,0))
boxplot(log2(expf+1))

###Normalize
bolFMat<-as.matrix(expf, nrow = nrow(expf), ncol = ncol(expf))
library(NormalyzerDE)
expf2<-performQuantileNormalization(bolFMat, noLogTransform = T)
rownames(expf2)<-rownames(expf)
colnames(expf2)<-colnames(expf)
expf2<-as.data.frame(expf2)
boxplot(log2(expf2+1))
write.table(expf2, file = "Expression_sinica_normalized.txt", quote=F)
```
