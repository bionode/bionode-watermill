'use strict'

// === WATERWHEEL ===
const waterwheel = require('./lib/waterwheel')
const { task, join, parallel } = waterwheel
const { File } = waterwheel.types
const { shell, shellPipe } = waterwheel.wrappers

// === PIPELINE ===
const fs = require('fs')
const ncbi = require('bionode-ncbi')
const request = require('request')

const THREADS = 4
const config = {
  sraAccession: '2492428',
  referenceURL: 'http://ftp.ncbi.nlm.nih.gov/genomes/all/GCA_000988525.2_ASM98852v2/GCA_000988525.2_ASM98852v2_genomic.fna.gz'
}

// task(props, cb)
// cb(props) will be called after Files in props are resolved with globby
// cb is the "task action creator" and should return a Stream, Promise, or callback

const samples = task({
  input: {
    db: 'sra',
    accession: config.sraAccession
  },
  output: '**/*.sra'
}, ({ input }) => ncbi.download(input.db, input.accession) )

const fastqDump = task({
  input: new File('**/*.sra'),
  output: [1, 2].map(n => new File(`*_${n}.fastq.gz`))
}, ({ input }) => shell(`fastq-dump --split-files --skip-technical --gzip ${input}`) )

const downloadReference = task({
  input: config.referenceURL,
  output: new File(config.referenceURL.split('/').pop())
}, ({ input, output }) => request(input).pipe(fs.createWriteStream(output.value)) )

const bwaIndex = task({
  input: new File('*_genomic.fna.gz'),
  output: ['amb', 'ann', 'bwt', 'pac', 'sa'].map(suffix => new File(`*_genomic.fna.gz.${suffix}`))
}, ({ input }) => shell(`bwa index ${input}`) )

const alignAndSort = task({
  input: [new File('*_genomic.fna.gz'), new File('*.fastq.gz')],
  output: new File('reads.bam')
}, ({ input }) => shellPipe([
  `bwa mem -t ${THREADS} ${input[0]} ${input[1]} ${input[2]}`,
  'samtools view -Sbh -',
  'samtools sort - -o reads.bam'
]) )

const samtoolsIndex = task({
  input: new File('*.bam'),
  output: new File('*.bam.bai')
}, ({ input }) => shell(`samtools index ${input}`) )

const decompressReference = task({
  input: new File('*_genomic.fna.gz'),
  output: new File('*_genomic.fna')
}, ({ input }) => shell(`bgzip -d ${input}`) )

const mpileupAndCall = task({
  input: [new File('*_genomic.fna'), new File('*.bam'), new File('*.bam.bai')],
  output: new File('variants.vcf')
}, ({ input }) => shellPipe([
  `samtools mpileup -uf ${input[0]} ${input[1]}`,
  'bcftools call -c - > variants.vcf'
]) )

const pipeline = parallel({
  taskLists: [[samples], [downloadReference, bwaIndex]],
  next: join(alignAndSort, samtoolsIndex, decompressReference, mpileupAndCall)
})
pipeline()
