'use strict'

// === WATERWHEEL ===
const {
  task,
  join,
  shell
} = require('../..')

// === MODULES ===
const fs = require('fs')

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
}, ({ params }) => {
  const { url } = params
  const outfile = url.split('/').pop()

  // essentially curl -O
  return request(url).pipe(fs.createWriteStream(outfile))
})

// To run one task by itself:
// getReference()
//   .on('task.finish', (task) => {
//     console.log('output: ', task.output)
//   })

/**
 * Indexes a reference genome with bwa.
 * @input {string} the reference genome ending in _genomic.fna.gz
 * @output {file[]} the annotation, index, etc. files produced by "bwa index"
 * @action {shell}
 */
const bwaIndex = task({
  input: '*_genomic.fna.gz',
  output: ['amb', 'ann', 'bwt', 'pac', 'sa'].map(suffix => `*_genomic.fna.gz.${suffix}`),
  name: 'bwa index *_genomic.fna.gz',
}, ({ input }) => shell(`bwa index ${input}`) )

// bwaIndex()
//   .on('task.finish', function(task) {
//     console.log('output: ', task.resolvedOutput)
//   })

// Can join two tasks as long as the input for task 2 comes from output of task
// 1
// join(getReference, bwaIndex)()
//   .on('destroy', function() {
//     console.log('output: ', JSON.stringify(this._output, null, 2))
//   })
//
//   OR
//   to take input from fs:
//  task1().on('destroy', () => {
//    task2().on('destroy', ...)
//  })


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
  name: `Download SRA ${config.sraAccession}`
// }, ({ params }) => ncbi.download(params.db, params.accession) )
}, ({ params }) => ncbi.download(params.db, params.accession).resume() )

// getSamples()
//   .on('task.finish', function(task) {
//     console.log('output: ', task.resolvedOutput)
//   })

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
}, ({ input }) => shell(`fastq-dump --split-files --skip-technical --gzip ${input}`) )

// fastqDump()
//   .on('task.finish', (task) => {
//     console.log('output: ', task.resolvedOutput)
//   })

// join(getSamples, fastqDump)()
//   .on('task.finish', function(task) {
//     console.log('output: ', task.resolvedOutput)
//   })

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
    }
  },
  output: ['reads.bam'], // NOTE forced genomic.fna into here
  name: 'bwa mem | samtools view | samtools sort'
}, ({ input }) => shell(`
bwa mem -t ${THREADS} ${input.reference} ${input.reads['1']} ${input.reads['2']} | \
samtools view -@ ${THREADS} -Sbh - | \
samtools sort -@ ${THREADS} - -o reads.bam > reads.bam
`))

join(getReference, bwaIndex, getSamples, fastqDump, alignAndSort)()
  .on('task.finish', function(task) {
    console.log('task: ', task)
  })

// alignAndSort()
//   .on('destroy', function() {
//     console.log('output: ', JSON.stringify(this._output, null, 2))
//   })

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
}, ({ input }) => shell(`samtools index ${input}`))


// join(alignAndSort, samtoolsIndex)()
//   .on('destroy', function() {
//     console.log('output: ', JSON.stringify(this._output, null, 2))
//   })

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
}, ({ input }) => shell(`bgzip -d ${input} --stdout > ${input.slice(0, -('.gz'.length))}`) )


// join(getReference, decompressReference, alignAndSort, samtoolsIndex)()
//   .on('destroy', function() {
//     console.log('output: ', JSON.stringify(this._output, null, 2))
//   })

/**
 * Multipileup bam and call variants.
 * @input.reference {file} the reference genome
 * @input.bam {file} the binary alignment map
 * @input._bai {file} the binary alignment index, not used directly but required
 * @output {file} variant calling format
 * @action {shell}
 */
const mpileupandcall = task({
  input: {
    reference: '*_genomic.fna',
    bam: '*.bam',
    _bai: '*.bam.bai'
  },
  output: 'variants.vcf',
  name: 'samtools mpileup | bcftools call',
  // skipPassed: true
}, ({ input }) => shell(`
samtools mpileup -uf ${input.reference} ${input.bam} | \
bcftools call -c - > variants.vcf
`))

// this fails b/c mpileupandcall needs reference *_genomic.fna and wont
// find it in the output of samtoolsindex
// join(decompressReference, alignAndSort, samtoolsIndex, mpileupandcall)()
//   .on('destroy', function() {
//     console.log('output: ', JSON.stringify(this._output, null, 2))
//   })

// === task orchestration ===

// const reads = join(getSamples, fastqDump)
// const reference = join(getReference, parallel(bwaIndex, decompressReference))
// const call = join(alignAndSort, samtoolsIndex, mpileupandcall)
//
// const pipeline = join(parallel(reads, reference), call)
//
// call()
//   .on('close', function() {
//     this.output()
//     console.log('Pipeline has finished')
//   })
