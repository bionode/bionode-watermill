# Comparisons

WIP

bionode-watermill is an API targeted towards developing data analysis
pipelines. It separates itself from other tools by putting **parameters ahead
of input/output files**. This enables modular rearrangement of pipelines, and
quick experimentation and comparison of tools and their parameters. Let's
examine how other tools put I/O ahead of parameters, then demonstrate how
bionode-watermill does the opposite.

Consider a simple variant calling pipeline:

1. Download reference genome
2. Index reference genome
3. Download reads
4. Extract reads into fastq
5. Align fastq to reference producing alignment in BAM format
6. Index bam
7. Decompress reference (only b/c mpileup cannot take *.fna.gz)
8. mpileup | call to produce variants in VCF format

As a makefile:
For those not versed in make:

```
target: prerequisites
    command
$@ is target, $< is prerequisites
```

TODO replace this with snakemake.

```make
READS=2492428
REFERENCE_URL=http://ftp.ncbi.nlm.nih.gov/genomes/all/GCA_000988525.2_ASM98852v2/GCA_000988525.2_ASM98852v2_genomic.fna.gz

# make and snakemake are "pull" style: declare final files, then declare
# tasks to make them, until you "go backwards" to the "source". In this
# pipeline, the sources is the READS and REFERENCE_URL variables.
# The final file this pipeline will produce is variants.vcf
all: variants.vcf

# Step 1
reference.genomic.fna.gz:
  curl -s ${REFERENCE_URL} -o $@

# Step 2
reference.genomic.fna.gz.bwt: reference.genomic.fna.gz
  bwa index $<

# Step 3
reads.sra:
  bionode-ncbi download sra ${READS} > /dev/null; \
  cp ${READS}/*.sra $@

# Step 4
fastq-dump.log: reads.sra
  fastq-dump --split-files --skip-technical --gzip $< > $@

# Step 5
reads.bam: reference.genomic.fna.gz reference.genomic.fna.bwt fastq-dump.log
  bwa mem reference.genomic.fna.gz reads_1.fastq.gz reads_2.fastq.gz | \
  samtools view -Sbh - | \
  samtools sort - -o reads.bam

# Step 6
reads.bam.bai: reads.bam
  samtools index $<

# Step 7
reference.genomic.fna: reference.genomic.fna.gz
  bgzip -d $< --stdout > $@

# Step 8
variants.vcf: reference.genomic.fna reads.bam
  samtools mpileup -uf reference.genomic.fna reads.bam | \
  bcftools call -c - > $@
```

You can notice that each task is described by the prerequisite files and the
target file to produce. Command parameters and the tools themselves are
nested internally. Different tools produce different files: switching tools
around would require modifying targets and prerequisites in the preceding
and following tasks.

Although Snakemake is much more elegant, and good use of its wildcards can
make for a succinct pipeline, the same issues still persist. Nextflow uses
a "dataflow variable" approach which is an improvement. Note also that
Nextflow uses a "push" style, where a task produces files by inserting them
into queues called "channels", and then downstream tasks take input from
those channels. A portion of the pipeline using Nextflow:

```groovy
#!/usr/bin/env
READS = 2492428
REFERENCE_URL = http://ftp.ncbi.nlm.nih.gov/genomes/all/GCA_000988525.2_ASM98852v2/GCA_000988525.2_ASM98852v2_genomic.fna.gz

process downloadReference {
    input: val referenceURL from REFERENCE_URL
    output: file 'reference.genomic.fna.gz' into referenceGenomeGz

    """
    appropriate/curl $referenceURL -o reference.genomic.fna.gz
    """
}

referenceGenomeGz.into { referenceGenomeGz1;
                         referenceGenomeGz2;
                         referenceGenomeGz3;
                         referenceGenomeGz4 }

process downloadSRA {
    input: val readsID from READS
    output: file '**.sra' into reads

    """
    bionode-ncbi download sra $readsID > tmp
    """
}
```
