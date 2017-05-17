'use strict'

// === WATERWHEEL ===
const {
  task,
  join,
  junction,
} = require('../../..')

// === MODULES ===
const fs = require('fs-extra')
const path = require('path')
const through = require('through2')

const request = require('request')
const ncbi = require('bionode-ncbi')

// === CONFIG ===
const THREADS = parseInt(process.env.WATERWHEEL_THREADS) || 4
const config = {
  name: 'Salmonella enterica',
  sraAccession: '2492428',
  referenceURL: 'http://ftp.ncbi.nlm.nih.gov/genomes/all/GCA_000988525.2_ASM98852v2/GCA_000988525.2_ASM98852v2_genomic.fna.gz'
}

// === TASKS ===

/**
 * Downloads the reference genome from NCBI.
 * @input {string} the reference url
 * @output {file} the reference genome in fna.gz
 * @action http request | writable
 */
const getReference = task({
  params: { url: config.referenceURL },
  input: null,
  output: '*_genomic.fna.gz',
  name: `Download reference genome for ${config.name}`
}, ({ params, dir }) => {
  const { url } = params
  const outfile = url.split('/').pop()

  // essentially curl -O
  return request(url).pipe(fs.createWriteStream(dir + '/' + outfile))
})


/**
 * Indexes a reference genome with bwa.
 * @input {string} the reference genome ending in _genomic.fna.gz
 * @output {file[]} the annotation, index, etc. files produced by "bwa index"
 * @action {shell}
 */
const bwaIndex = task({
  input: '*_genomic.fna.gz',
  output: ['amb', 'ann', 'bwt', 'pac', 'sa'].map(suffix => `*_genomic.fna.gz.${suffix}`),
  name: 'bwa index *_genomic.fna.gz'
}, ({ input }) => `bwa index ${input}` )


/**
 * Downloads the reads samples in SRA format.
 * @input.db {string} the NCBI database
 * @input.accession {string} the accession ID
 * @output {file} the reads
 * @action {readable} the stream returned by ncbi.download
 */
const getSamples = task({
  params: {
    db: 'sra',
    accession: config.sraAccession
  },
  input: null,
  output: '**/*.sra',
  dir: process.cwd(),
  name: `Download SRA ${config.sraAccession}`
}, ({ params }) => ncbi.download(params.db, params.accession).resume() )


/**
 * Extracts SRA into fastq using fastq-dump.
 * @input {file} the reads
 * @output {file[]} the written ends
 * @action {shell}
 */
const fastqDump = task({
  input: '**/*.sra',
  output: [1, 2].map(n => (`*_${n}.fastq.gz`)),
  name: 'fastq-dump **/*.sra'
}, ({ input }) => `fastq-dump --split-files --skip-technical --gzip ${input}` )


/**
 * Align reads and sort.
 * @input.reference {file} reference genome
 * @input.reads.1 {file} reads end 1
 * @input.reads.2 {file} reads end 2
 * @output {file} bam
 * @action {shell}
 */
const alignAndSort = task({
  input: {
    reference: '*_genomic.fna.gz',
    reads: {
      1: '*_1.fastq.gz',
      2: '*_2.fastq.gz'
    },
    indexFiles: ['amb', 'ann', 'bwt', 'pac', 'sa'].map(suffix => `*_genomic.fna.gz.${suffix}`)
  },
  output: ['reads.bam'], // NOTE forced genomic.fna into here
  name: 'bwa mem | samtools view | samtools sort'
}, ({ input }) => `
bwa mem -t ${THREADS} ${input.reference} ${input.reads['1']} ${input.reads['2']} | \
samtools view -@ ${THREADS} -Sbh - | \
samtools sort -@ ${THREADS} - -o reads.bam > reads.bam
`)


/**
 * Index bam.
 * @input {file} the bam
 * @output {file} the binary alignment index
 * @action {shell}
 */
const samtoolsIndex = task({
  input: 'reads.bam',
  output: '*.bam.bai',
  name: 'samtools index'
}, ({ input }) => `samtools index ${input}`)


/**
 * Decompress a reference genome.
 * @input {file} the reference genome
 * @output {file} the decompressed reference genome
 * @action {shell}
 */
const decompressReference = task({
  input: '*_genomic.fna.gz',
  output: '*_genomic.fna',
  name: 'Decompress reference'
}, ({ input }) => `bgzip -d ${input} --stdout > ${input.slice(0, -('.gz'.length))}` )


/**
 * Multipileup bam and call variants.
 * @input.reference {file} the reference genome
 * @input.bam {file} the binary alignment map
 * @input._bai {file} the binary alignment index, not used directly but required
 * @output {file} variant calling format
 * @action {shell}
 */
const mpileupAndCall = task({
  input: {
    reference: '*_genomic.fna',
    bam: '*.bam',
    _bai: '*.bam.bai'
  },
  output: 'variants.vcf',
  name: 'samtools mpileup | bcftools call',
}, ({ input }) => `
samtools mpileup -uf ${input.reference} ${input.bam} | \
bcftools call -c - > variants.vcf
`)

// === PIPELINE ===

const pipeline = join(
  junction(
    join(getReference, bwaIndex),
    join(getSamples, fastqDump)
  ),
  decompressReference, // only b/c mpileup did not like fna.gz
  join(alignAndSort, samtoolsIndex, mpileupAndCall)
)

pipeline().then(results => console.log('PIPELINE RESULTS: ', results))

