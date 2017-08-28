'use strict'

/* DISCLAIMER: This is a very experimental pipeline and currently it is not
* solved. If you have an idea on how to solve this problem please make a PR.
* The idea with this pipeline is for the two-mappers pipeline.js to be able
* to process several samples.
*/

// === WATERMILL ===
const {
  task,
  join,
  junction,
  fork
} = require('../../..')

// === MODULES ===

const fs = require('fs')
const path = require('path')

const request = require('request')

// === CONFIG ===

// specifies the number of threads to use by mappers
const THREADS = parseInt(process.env.WATERMILL_THREADS) || 2

const config = {
  name: 'Streptococcus pneumoniae',
  sraAccession: ['ERR045788', 'ERR016633'],
  referenceURL: 'http://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/007/045/GCF_000007045.1_ASM704v1/GCF_000007045.1_ASM704v1_genomic.fna.gz'
}

// === TASKS ===

// first lets get the reference genome for our mapping
const getReference = task({
  params: { url: config.referenceURL },
  output: '*_genomic.fna.gz',
  dir: process.cwd(),   // this sets directory for outputs!
  name: `Download reference genome for ${config.name}`
}, ({ params, dir }) => {
  const { url } = params
  const outfile = url.split('/').pop()

  // essentially curl -O
  return request(url).pipe(fs.createWriteStream(outfile))

})

//then get samples to work with
const getSamples = (sraAccession) => task({
  params: {
    db: 'sra',
    accession: sraAccession
  },
  output: '**/*.sra',
  dir: process.cwd(), // Set dir to resolve input/output from
  name: `Download SRA ${sraAccession}`
}, ({ params }) => `bionode-ncbi download ${params.db} ${params.accession}`
)

// extract the samples from fastq.gz
const fastqDump = task({
  input: '**/*.sra',
  output: [1, 2].map(n => `*_${n}.fastq.gz`),
  name: 'fastq-dump **/*.sra'
}, ({ input }) => `fastq-dump --split-files --skip-technical --gzip ${input}`
)

// first lets uncompress the gz
const gunzipIt = task({
    input: '*_genomic.fna.gz',
    output: '*.fna',
    name: 'gunzipIt task'
  }, ({ input }) => `gunzip -c ${input} > ${input.split('.').slice(0,2).join('.')}.fna`
)

// then index using first bwa ...
const indexReferenceBwa = task({
  input: '*.fna',
  output: {
    indexFile: ['amb', 'ann', 'bwt', 'pac', 'sa'].map(suffix =>
      `bwa_index.${suffix}`),
    //reference: 'bwa_index.fa' //names must match for bwa - between reference
    // and index files
  },
  //params: { output: 'bwa_index.fa' },
  name: 'bwa index bwa_index.fna -p bwa_index'
}, ({ input }) => `bwa index ${input} -p bwa_index`)

// and bowtie2

const indexReferenceBowtie2 = task({
    input: '*.fna',
    output: ['1.bt2', '2.bt2', '3.bt2', '4.bt2', 'rev.1.bt2',
      'rev.2.bt2'].map(suffix => `bowtie_index.${suffix}`),
    params: { output: 'bowtie_index' },
    name: 'bowtie2-build -q uncompressed.fna bowtie_index'
  }, ({ params, input }) => `bowtie2-build -q ${input} ${params.output}`
  /* for bowtie we had to uncompress the .fna.gz file first before building
   the reference */
)

// now use mappers with bwa

const bwaMapper = task({
    input: {
      //reference: '*.fa',
      reads:[1, 2].map(n => `*_${n}.fastq.gz`),
      indexFiles: ['amb', 'ann', 'bwt', 'pac', 'sa'].map(suffix =>
        `bwa_index.${suffix}`) //pass index files to bwa mem
    },
    output: '*.sam',
    params: { output: 'bwa_output.sam' },
    name: 'Mapping with bwa...'
  }, ({ input, params }) => `bwa mem -t ${THREADS} bwa_index ${input.reads[0]} ${input.reads[1]} > ${params.output}`
)

// and with bowtie2

const bowtieMapper = task({
  input: {
    reference: '*_genomic.fna.gz',
    reads:[1, 2].map(n => `*_${n}.fastq.gz`),
    indexFiles: ['1.bt2', '2.bt2', '3.bt2', '4.bt2', 'rev.1.bt2',
      'rev.2.bt2'].map(suffix => `bowtie_index.${suffix}`) //pass index files to
    // bowtie2
  },
  output: '*.sam',
  params: { output: 'bowtie2_output.sam' },
  name: 'Mapping with bowtie2...'
}, ({ input, params }) => `bowtie2 -p ${THREADS} -x bowtie_index -1 ${input.reads[0]} -2 ${input.reads[1]} -S ${params.output}`
)

// === PIPELINE ===

// first gets reference
getReference().then(results => {
//  console.log("results:", results.resolvedOutput)
  const pipeline = (sraAccession) => join(
    getSamples(sraAccession),
    fastqDump,
    gunzipIt,
    fork(
      join(indexReferenceBwa, bwaMapper),
      join(indexReferenceBowtie2, bowtieMapper)
    )
  )
// then fetches the samples and executes the remaining pipeline
  for (const sra of config.sraAccession) {
    //console.log("sample:", sra, results.resolvedOutput)
    console.log("cwd_check", process.cwd())
    const pipelineMaster = pipeline(sra)
    pipelineMaster().then(results => console.log("Results: ", results))
  }
})
