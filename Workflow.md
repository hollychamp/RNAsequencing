Convert .gff files provided by microbesNG to .gtf files so they can be used with STAR

<details>
<summary> Note </summary>
This .gtf file does not carry over uniprotKB annotations but they will be merged in later steps in R
</details>

```
gffread genome1.gff -T -o genome1.gtf

```

Use STAR to map bacterial transcripts to bacterial genomes 

Build genome index

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
