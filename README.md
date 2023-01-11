# The Impact of Spatial Correlation on Methylation Entropy with Application to Mouse Brain Methylome
This is a source code for analysis in "The Impact of Spatial Correlation on Methylation Entropy with Application to Mouse Brain Methylome" paper.

## Requirements
* SRA toolkit (>= 3.0.1)
* TrimGalore (>= 0.6.7)
* Bismark with bowtie2 (>= 2.4.4)
* Samtools (>= 1.4)

## 1. Download raw FastQ files from GEO database
Take [GSE47966](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE47966) as example, to get the FastQ files for the 1st sample :

1) Go to the webpage for that sample by clicking [GSM1163695](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM1163695).
2) Click “SRA Run Selector” at the bottom to get related SRA information (what we needed is the SRA accession for each data file)
3) Identify associated SRA accessions (eg. SRR901379, SRR901380 etc), then run fastq-dump for each of them in the terminal.
```
fastq-dump —split-3 SRR901379
```

## 2. Trim bad quality bases using Trim Galore
```
trim_galore SRR901379.fastq
```

## 3. Generate the index of reference genome
[Bismark](https://www.bioinformatics.babraham.ac.uk/projects/bismark/) is needed, which depends on other tools including Bowtie2 and samtools. 
1) Download reference genome of FastA format. For the latest release of human genome (hg38), it can be found from [this link](http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/). To download the file, type the following command in the terminal:
```
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
```
2) Uncompress hg38.fa.gz to get hg38.fa
```
gunzip hg38.fa.gz
```
3) Put hg38.fa to a folder hg38Idx, this folder will contain indexed genome
```
mkdir hg38Idx
mv hg38.fa hg38Idx
```
4) Generate genome index
```
bismark_genome_preparation hg38Idx
```
All output file will be generated in the folder hg38Idx.

## 4. Align reads to reference genome
1) Align reads to reference genome
```
bismark —parallel 4 -p 4 hg38Idx SRR901379_trimmed.fq
```
* The output file will be: SRR901379_trimmed_bismark_bt2.bam

2) Discard PCR duplicate
```
deduplicate_bismark -s -o SRR901379_trimmed_bismark SRR901379_trimmed_bismark_bt2.bam 
```
* The output file will be: hg38Idx SRR901379_trimmed_bismark.deduplicated.bam

3) Perform methylation calling
```
bismark_methylation_extractor -s —comprehensive —merge_non_CpG —parallel 6 —genome_folder ./hg38Idx SRR901379_trimmed_bismark.deduplicated.bam
```
The output file with CpG methylation pattern is: CpG_context_XXX_bismark_bt2.txt contain methylation calling for CpGs. Such information will be merged to get patterns for 4CpG segments.

## 5. Extract methylation pattern for 4CpG segments
Clone the repository or download source code files.

1) Get the coordinates for all CpG site from genome:
```
perl get_string_position_from_genome.pl hg38.fa CG >hg38.CpG.pos.txt
```

2) Get detailed methylation pattern for each CpG site:
```
perl merge_CpG_methyl_info.pl ./hg38.CpG.pos.txt CpG_context_SRR901379_trimmed_bismark.deduplicated.txt >SRR901379_trimmed_bismark.CpG.inf.txt
```

3) Get detailed methylation pattern for each 4CG segment:
```
perl extract_4CG_seed_methyl_info.pl ~/db/CpG/hg38.CpG.pos.txt SRR901379_trimmed_bismark.CpG.inf.txt >SRR901379_trimmed_bismark.CpG.4CG.txt
```

4) Convert the 4CG methyl file to the desired format
```
perl convert_4CG_seed_format.pl SRR901379_trimmed_bismark.CpG.4CG.txt >SRR901379_trimmed_bismark.CpG.4CG.cnv.txt
```

5) Further filtering based on read depth etc. For example, to keep 4CG segments with depth >=10:
```
cat SRR901379_trimmed_bismark.CpG.4CG.cnv.txt | awk ‘$3>=10’ > SRR901379_trimmed_bismark.CpG.4CG.cnv.10X.txt
```

* The final output file is of format with the header line shows the content for each column. Each line is the information for one 4CG segment. Importantly, the last column is the methylation pattern for the 4CG segment, including the counts for each pattern. For example, 0000:7;1110:2;1111:2 means there are three patterns (0000 means all 4 CpG are unmethylated, 1110 means CpG 1-3 are methylated and CpG4 unmethylated, 1111 means all 4 CpG are methylated), with the number of reads supporting each pattern are 7, 2 and 2, respectively.

## R scripts for Simulation 1, 2, and Real data analysis in the manuscript
* Simulation1.R
* Simulation2.R 
* RealDemonstration.R 

## Contact
If you have any questions or problems, please contact to joungmin AT vt.edu.
