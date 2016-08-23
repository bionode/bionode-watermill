'use strict'

// === WATERMILL ===
// TODO npm publish and require bionode-watermill
const {
  task,
  join,
  junction,
  fork
} = require('../..')

// === MODULES ===
const fs = require('fs')
const path = require('path')

const request = require('request')
const ncbi = require('bionode-ncbi')

// === CONFIG ===
const THREADS = parseInt(process.env.WATERMILL_THREADS) || 4
const MEMORYGB = parseInt(process.env.WATERMILL_MEMORY) || 4
const TMPDIR = path.resolve(__dirname, 'temp') // Assume this already exists
const config = {
  name: 'Salmonella enterica',
  sraAccession: '2492428',
  referenceURL: 'http://ftp.ncbi.nlm.nih.gov/genomes/all/GCA_000988525.2_ASM98852v2/GCA_000988525.2_ASM98852v2_genomic.fna.gz'
}

const KMERSIZE = 20
const MINCOVERAGE = 5
const PLOTXMAX = 60
const PLOTYMAX = 1200000

// === TASKS ===

/**
 * Downloads the reference genome from NCBI.
 * @input {string} the reference url
 * @output {file} the reference genome in fna.gz
 * @action http request | writable
 */
const getReference = task({
  params: { url: config.referenceURL },
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
  output: '**/*.sra',
  dir: process.cwd(), // Set dir to resolve input/output from
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
  output: [1, 2].map(n => `*_${n}.fastq.gz`),
  name: 'fastq-dump **/*.sra'
}, ({ input }) => `fastq-dump --split-files --skip-technical --gzip ${input}` )


/**
 * Trimming with trimmomatic
 */
const trim = task({
  input: [1, 2].map(n => `*_${n}.fastq.gz`),
  output: [1, 2].map(n => `*_${n}.trim.pe.fastq.gz`),
  params: {
    trimmomaticExecutable: path.resolve(__dirname, 'bin', 'trimmomatic-0.36.jar'),
    adapters: path.resolve(__dirname, 'adapters')
  },
  name: 'Trimming with trimmomatic'
}, ({ input, params }) => `java -jar ${params.trimmomaticExecutable} PE -phred33 \
${input[0]} ${input[1]} \
reads_1.trim.pe.fastq.gz reads_1.trim.se.fastq.gz \
reads_2.trim.pe.fastq.gz reads_2.trim.se.fastq.gz \
ILLUMINACLIP:${params.adapters}/TruSeq3-PE.fa:2:30:10 \
LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36`)


/**
 * Merge trim ends with seqtk merge.
 */
const mergeTrimEnds = task({
  input: [1, 2].map(n => `*_${n}.trim.pe.fastq.gz`),
  output: 'reads.trim.pe.fastq.gz',
  name: 'Merge trim ends'
}, ({ input }) => `seqtk mergepe ${input[0]} ${input[1]} | gzip > reads.trim.pe.fastq.gz`)


/**
 * Fitering with KMC.
 */
const filterKMC = task({
  input: 'reads.trim.pe.fastq.gz',
  output: 'reads.trim.pe.filtered.fastq.gz',
  params: { kmcFile: 'reads.trim.pe.kmc' },
  name: 'Filtering with KMC'
}, ({ input, params }) => `
kmc -k${KMERSIZE} -m${MEMORYGB} -t${THREADS} ${input} ${params.kmcFile} ${TMPDIR} && \
kmc_tools filter ${params.kmcFile} -cx${MINCOVERAGE} ${input} -ci0 -cx0 reads.trim.pe.filtered.fastq.gz
`)


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
  trim, mergeTrimEnds, filterKMC
  // decompressReference, // only b/c mpileup did not like fna.gz
  // join(alignAndSort, samtoolsIndex, mpileupAndCall)
)

pipeline().then(results => console.log('PIPELINE RESULTS: ', results))

