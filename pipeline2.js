'use strict'

const fs = require('fs')
const request = require('request')
const ncbi = require('bionode-ncbi')

const Task = require('./lib/Task.js')
const waterwheel = require('./lib/waterwheel')
const { shell, shellPipe, Process } = waterwheel.wrappers
const Join = require('./lib/Join.js')

const THREADS = 4
const config = {
  name: 'Salmonella enterica',
  sraAccession: '2492428',
  referenceURL: 'http://ftp.ncbi.nlm.nih.gov/genomes/all/GCA_000988525.2_ASM98852v2/GCA_000988525.2_ASM98852v2_genomic.fna.gz'
}

// Reference

const downloadReference = Task({
  input: { value: config.referenceURL },
  output: { file: config.referenceURL.split('/').pop() },
  name: `Download reference genome for ${config.name}` 
}, ({ input }) => request(input).pipe(fs.createWriteStream(input.split('/').pop())) )

const bwaIndex = Task({
  input: { file: '*_genomic.fna.gz' },
  output: ['amb', 'ann', 'bwt', 'pac', 'sa'].map(suffix => ({ file: `*_genomic.fna.gz.${suffix}` })),
  name: 'bwa index *_genomic.fna.gz'
}, ({ input }) => new Process(`bwa index ${input}`) )


const pipeline = Join(downloadReference, bwaIndex)

pipeline()
  .on('task.done', (output) => {
    console.log('Join emitted a task.done')
    console.log(output)
  })

// Samples

// const samples = Task({
//   input: {
//     db: 'sra',
//     accession: config.sraAccession
//   },
//   output: new File('*|)}>#*.sra'),
//   name: `Download SRA ${config.sraAccession}`
// }, ({ input }) => ncbi.download(input.db, input.accession) )
//
// samples()
