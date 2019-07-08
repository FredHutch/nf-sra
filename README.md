# Nextflow - SRA
Nextflow pipeline to download FASTQ data from the NCBI SRA repository

### Executing

Run with the command:

```
nextflow run FredHutch/nf-sra [OPTIONS]
```


 ### Options

 To download a single SRA accession, use --sra <ACCESSION>.

 To download an entire BioProject, use --bioproject <ACCESSION>.

 In either case, use --output_folder to specify the location for the downloaded files.
 All files will be named for the SRA accession: <ACCESSION>.fastq.gz.
 Paired-end files will be interleaved.
 