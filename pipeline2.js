'use strict'

const fs = require('fs')
const request = require('request')

const Task = require('./lib/Task.js')
const waterwheel = require('./lib/waterwheel')
const { File } = waterwheel.types
const { shell, shellPipe } = waterwheel.wrappers
const Join = require('./lib/Join.js')

const THREADS = 4
const config = {
  name: 'Salmonella enterica',
  sraAccession: '2492428',
  referenceURL: 'http://ftp.ncbi.nlm.nih.gov/genomes/all/GCA_000988525.2_ASM98852v2/GCA_000988525.2_ASM98852v2_genomic.fna.gz'
}

const downloadReference = Task({
  input: config.referenceURL,
  output: new File(config.referenceURL.split('/').pop()),
  name: `Download reference genome for ${config.name}` 
}, ({ input, output }) => request(input).pipe(fs.createWriteStream(output.value)) )

const bwaIndex = Task({
  input: new File('*_genomic.fna.gz'),
  output: ['amb', 'ann', 'bwt', 'pac', 'sa'].map(suffix => new File(`*_genomic.fna.gz.${suffix}`)),
  name: 'bwa index *_genomic.fna.gz'
}, ({ input }) => shell(`bwa index ${input}`) )

// downloadReference()
//   .on('task.done', (output) => console.log(output))

// bwaIndex()

const pipeline = Join(downloadReference, bwaIndex)
