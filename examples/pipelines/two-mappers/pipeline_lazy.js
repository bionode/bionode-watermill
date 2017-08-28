'use strict'

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
const ncbi = require('bionode-ncbi')

// === CONFIG ===

// specifies the number of threads to use by mappers
const THREADS = parseInt(process.env.WATERMILL_THREADS) || 2

const config = {
  name: 'Streptococcus pneumoniae',
  sraAccession: 'ERR045788',
  referenceURL: 'http://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/007/045/' +
  'GCF_000007045.1_ASM704v1/GCF_000007045.1_ASM704v1_genomic.fna.gz'
}

// === TASKS ===

// first lets get the reference genome for our mapping
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

// then get samples to work with
const getSamples = task({
    params: {
      db: 'sra',
      accession: config.sraAccession
    },
    output: '**/*.sra',
    dir: process.cwd(), // Set dir to resolve input/output from
    name: `Download SRA ${config.sraAccession}`
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
  output: '*.fa',
  params: { output: 'uncompressed.fa' }
}, ({ params, input}) => `gunzip -c ${input} > ${params.output}`
)

// then index using first bwa ...
const indexReferenceBwa = task({
  input: '*.fa',
  output: {
    indexFile: ['amb', 'ann', 'bwt', 'pac', 'sa'].map(suffix =>
      `bwa_index.${suffix}`),
    //reference: 'bwa_index.fa' //names must match for bwa - between reference
    // and index files
  },
  //params: { output: 'bwa_index.fa' },
  name: 'bwa index bwa_index.fa -p bwa_index'
}, ({ input }) => `bwa index ${input} -p bwa_index`)

// and bowtie2

const indexReferenceBowtie2 = task({
    input: '*.fa',
    output: ['1.bt2', '2.bt2', '3.bt2', '4.bt2', 'rev.1.bt2',
      'rev.2.bt2'].map(suffix => `bowtie_index.${suffix}`),
    params: { output: 'bowtie_index' },
    name: 'bowtie2-build -q uncompressed.fa bowtie_index'
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

// SAMTOOLS tasks

const samtoolsFaidx = task({
  input: '*_genomic.fna.gz',
  output: '*.fasta.fai',
  name: 'generating .fai for reference'
}, ({ input }) => `gunzip -c ${input} > reference.fasta && samtools faidx reference.fasta`
)

const samtoolsView = task({
  input: {
    samfile: '*.sam',
    faidxfile: '*.fasta.fai'
  },
  output: '*.bam',
  name: 'Samtools View'
}, ({ input }) => {
  const outputfile = input.samfile.slice(0,-3) + 'bam'  /* use same name as
   input but replace the final extension to the string to bam. In this case
    if we hardcode the name, it will render an error since two outputs will
     match the same name */
  return `samtools view -b -S -t ${input.faidxfile} -@ ${THREADS} -o ${outputfile} ${input.samfile}`
})

const samtoolsSort = task({
  input: '*.bam',
  output: '*.sorted.bam'
}, ({ input }) => {
  const outputsorted = input.slice(0,-3) + 'sorted.bam' /* use same name as
   input but replace the final extension to the string to sorted.bam. */
  return `samtools sort -@ ${THREADS} -o ${outputsorted} ${input}`
})

const samtoolsIndex = task({
  input: '*.sorted.bam',
  output: '*.bai'
}, ({ input }) => `samtools index ${input}` /* this will generate a .bai
 file based on the input file name and thus there is no need to handle the
  file extension for the next task */
)

const samtoolsDepth = task({
  input: {
    indexfile: '*.bai', /* index is needed although it is not an input to
     the command of samtools depth */
    sortedfile: '*.sorted.bam'
  },
  output: '*.txt'
}, ({ input }) => {
  const depthoutput = input.sortedfile.slice(0,-11) + '_depth.txt'
  return `samtools depth ${input.sortedfile} > ${depthoutput}`
})


// === PIPELINE ===

const pipeline = join(
  junction(
    getReference,
    join(getSamples,fastqDump)
  ),
  samtoolsFaidx, gunzipIt, /* since this will be common to both mappers, there
   is no
   need to be executed after fork duplicating the effort. */
  fork(
    join(indexReferenceBwa, bwaMapper),
    join(indexReferenceBowtie2, bowtieMapper)
  ),
  /* here the pipeline will break into two distinct branches, one for each
   mapping approach used. */
  samtoolsView, samtoolsSort, samtoolsIndex, samtoolsDepth
)

// actual run pipelines and return results
pipeline().then(results => console.log('PIPELINE RESULTS: ', results))
