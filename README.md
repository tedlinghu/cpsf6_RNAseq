# cpsf6_RNAseq

<b>Daphne Cornish</b>, other authors, <b>Ted Ling Hu</b>,<b>Judd F. Hultquist</b>

<hr>

This repository contains the scripts needed to generate the figures and analysis as reported in Cornish et al (Unpublished). The script may need to be adapted to the local environment. Due to IRB constraints we are unable to share clinical data used to generate this data.


# Highlights
<hr>
<ul>
  <li>HIGHLIGHT1</li>
  <li>HIGHLIGHT2</li>
  <li>HIGHLIGHT3</li>
</ul>

# Summary
<hr>
SUMMARYSUMMARYSUMMARYSUMMARYSUMMARYSUMMARYSUMMARYSUMMARYSUMMARYSUMMARYSUMMARYSUMMARYSUMMARYSUMMARYSUMMARYSUMMARYSUMMARYSUMMARYSUMMARYSUMMARYSUMMARY


# Dependencies
<hr>
R
<ul>
  <li> DESeq2 </li>
  <li> apeglm </li>
  <li> dplyr </li>
  <li> EnhancedVolcano </li>
  <li> clusterProfiler </li>
  <li> org.Hs.eg.db </li>
  <li> ReactomePA </li>
  <li> pheatmap </li>
  <li> AnnotationDbi </li>

</ul>


# Bioinformatic analyses (based off <a href="https://www.nature.com/articles/nprot.2016.095#Sec11">Pertea et al., 2016</a>)


### Trimming

trimmomatic SE -phred33 -threads 32 ${infile} ${base}_trimmed.fastq.gz ILLUMINACLIP:TruSeq3-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

### Alignment

hisat2 --dta -p 8 -x genome -U "${extracted_name}_R1_001_trimmed.fastq.gz,${extracted_name}_R2_001_trimmed.fastq.gz" -S "${extracted_name}.sam

### Assembly and Merge

stringtie -e -B -G genome.gtf -o "${base}_stringtie_output/${base}.gtf" -l "${base}_stringtie_output/${base}" "$infile"
stringtie --merge -p 8 -G genome.gtf -o merged.gtf mergelist.txt
stringtie -e -p 8 -B -G merged.gtf -o ${base}_stringtie_output_merged/${base}.gtf -l ${base}_stringtie_output_merged/${base} ${infile}

### Preparation for DESeq2

python /software/stringtie/2.1.3/prepDE.py -i sample_IDs.txt

