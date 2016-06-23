'use strict'

// === WATERWHEEL ===
const waterwheel = require('./lib/waterwheel')
const { task, join, parallel } = waterwheel
const { File } = waterwheel.types
const { shell, shellPipe } = waterwheel.wrappers

// === UTILS ===
const last = (str, sep) => {
  const splitted = str.split(sep)
  return splitted[splitted.length - 1]
}

// === PIPELINE ===
const fs = require('fs')
const ncbi = require('bionode-ncbi')
const request = require('request')

const THREADS = 4
const config = {
  sraAccession: '2492428',
  referenceURL: 'http://ftp.ncbi.nlm.nih.gov/genomes/all/GCA_000988525.2_ASM98852v2/GCA_000988525.2_ASM98852v2_genomic.fna.gz'
}

// const samples = task(...)
// const samples = () => task(...)
function samples(next) {
  // task(props, cb)(next)
  return task(
    // these params are then made available to the cb
    {
      input: {
        db: 'sra',
        accession: config.sraAccession
      },
      // will be globby'ed after cb completes
      output: '**/*.sra'
    },
    // the cb for this task
    ( ({ input }) => ncbi.download(input.db, input.accession) )
  )(next)
}

function fastqDump(next) {
  return task(
  {
    input: new File('**/*.sra'),
    output: [1, 2].map(n => new File(`*_${n}.fastq.gz`))
  },
  ({ input }) => shell(`fastq-dump --split-files --skip-technical --gzip ${input}`)
)(next)
}

function downloadReference(next) {
  return task(
  {
    input: config.referenceURL,
    output: new File(last(config.referenceURL, '/'))
  },
  ({ input, output }) => request(input).pipe(fs.createWriteStream(output.value))
  )(next)
}

function bwaIndex(next) {
  return task(
  {
    input: new File('*_genomic.fna.gz'),
    output: ['amb', 'ann', 'bwt', 'pac', 'sa'].map(suffix => new File(`*_genomic.fna.gz.${suffix}`))
  },
  ({ input }) => shell(`bwa index ${input}`)
  )(next)
}

function alignAndSort(next) {
  return task(
  {
    input: [new File('*_genomic.fna.gz'), new File('*.fastq.gz')],
    output: new File('reads.bam')
  },
  ({ input }) => shellPipe([
    `bwa mem -t ${THREADS} ${input[0]} ${input[1]} ${input[2]}`,
    'samtools view -Sbh -',
    'samtools sort - -o reads.bam'
  ])
  )(next)
}

function samtoolsIndex(next) {
  return task(
  {
    input: new File('*.bam'),
    output: new File('*.bam.bai')
  },
  ({ input }) => shell(`samtools index ${input}`)
  )(next)
}

function decompressReference(next) {
  return task(
  {
    input: new File('*_genomic.fna.gz'),
    output: new File('*_genomic.fna')
  },
  ({ input }) => shell(`bgzip -d ${input}`)
  )(next)
}

function mpileupAndCall(next) {
  return task(
  {
    input: [new File('*_genomic.fna'), new File('*.bam'), new File('*.bam.bai')],
    output: new File('variants.vcf')
  },
  ({ input }) => shellPipe([
    `samtools mpileup -uf ${input[0]} ${input[1]}`,
    'bcftools call -c - > variants.vcf'
  ])
  )(next)
}

// function after(next) {
//   return task(
//   {
//     input: null,
//     output: null
//   },
//   () => shell('echo AFTER')
//   )(next)
// }

const pipeline = parallel({
  taskLists: [[samples], [downloadReference, bwaIndex]],
  // next: alignAndSort
  next: join(alignAndSort, samtoolsIndex, decompressReference, mpileupAndCall)
})
pipeline()
