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

const request = require('request')

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

//then get samples to work with
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

const samToolsView = task({
  input: '*.sam',
  output: '*.bam'
}, ({ input }) => `samtools view -Sb  ${input}  >  ${input.split('.')[0]}.bam`
)

const varCaller1 = task({
  input: {
    reference: '*.fna',
    bamFile: '*.bam'
  },
  output: '*1.vcf'
}, ({ input }) => `samtools mpileup -uf ${input.reference} ${input.bam} | bcftools call -c - > variants1.vcf`
)

const varCaller2 = task({
    input: {
      reference: '*.fna',
      bamFile: '*.bam'
    },
    output: '*2.vcf'
  }, ({ input }) => `samtools mpileup -uf ${input.reference} ${input.bam} | bcftools call -c - > variants2.vcf`
)

const varCaller3 = task({
    input: {
      reference: '*.fna',
      bamFile: '*.bam'
    },
    output: '*3.vcf'
  }, ({ input }) => `samtools mpileup -uf ${input.reference} ${input.bam} | bcftools call -c - > variants3.vcf`
)

const someOtherTask = task({
  input: {
    varCall1: '*1.vcf',
    varCall2: '*2.vcf',
    varCall3: '*3.vcf'
  },
  output: '.txt'
}, ({ input }) => `echo "${input.varCall1}\n${input.varCall2}\n${input.varCall3}" > listVarCalls.txt`)

// === PIPELINE ===

const pipeline = join(
  junction(
    getReference,
    join(getSamples,fastqDump)
  ),
  gunzipIt,
  indexReferenceBowtie2,
  bowtieMapper,
  samToolsView,
  junction(varCaller1, varCaller2, varCaller3),
  someOtherTask
)

// actual run pipelines and return results
pipeline().then(results => console.log('PIPELINE RESULTS: ', results))
