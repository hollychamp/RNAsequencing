**Step 1: Preparing .gtf files (optional)**

Convert .gff files provided by microbesNG to .gtf files so they can be used with STAR

<details>
<summary> Note </summary>
This .gtf file does not carry over uniprotKB annotations but they will be merged in later steps in R
</details>

<!-- This .gtf file does not carry over uniprotKB annotations but they will be merged in later steps in R -->

```
gffread genome1.gff -T -o genome1.gtf

```
**Step 2: Use STAR to map bacterial transcripts to bacterial genomes**

Build a genome index

```
STAR \
  --runMode genomeGenerate \
  --genomeDir APLE2372_index \
  --genomeFastaFiles genome1.fasta \
  --sjdbGTFfile genome1.gtf \
  --sjdbOverhang 6328143 \
  --runThreadN 40 \
  --genomeSAindexNbases 10

```

Map each set of transcripts to the genome 

```
STAR --runThreadN 40 \
  --genomeDir genome1_index \
  --readFilesIn transcript1_forward.fq.gz transcript1_reverse.fq.gz \
  --readFilesCommand gunzip \
  -c \
  --alignIntronMax 1 \
  --outFileNamePrefix ./genome1/genome1T1 \
  --quantMode GeneCounts \

```
