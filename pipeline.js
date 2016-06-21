'use strict'

// === WATERWHEEL ===
const { spawn } = require('child_process')
const globby = require('globby')

const task = (props, cb, next) => {
  const { output } = props

  const stream = cb(props)

  stream.on('end', () => globby(output).then(next) )

  return stream
}


// === PIPELINE ===
const ncbi = require('bionode-ncbi')

const config = {
  sraAccession: '2492428'
}

const samples = task(
  {
    input: {
      db: 'sra',
      accession: config.sraAccession
    },
    output: '**/*.sra'
  },
  ( ({ input }) => ncbi.download(input.db, input.accession) ),
  dump
  // ( (output) => dump(`fastq-dump --split-files --skip-technical --gzip ${output}`) )
)

// function dumps() {
//   return task({
//     input: '*.sra',
//     output: [1, 2].map(n => `*_${n}.fastq.gz`)
//   }, )
// }

function dump(data, cmd) {
  console.log(cmd)
  console.log(data[0])

  const process = spawn('fastq-dump', ['--split-files', '--skip-technical', '--gzip', data[0]])

  process.stdout.on('data', (data) => {
    console.log(`stdout: ${data}`)
  })

  process.stderr.on('data', (data) => {
    console.log(`stderr: ${data}`)
  })

  process.on('close', (code) => {
    console.log(`child process exited with code ${code}`)
  })
}

samples
  .on('data', console.log)
  .on('end', () => console.log('Finished SRA download'))

// ncbi.download('sra', config.sraAccession)
//   .on('data', console.log)
// samples
//   .on('end', path => console.log(`File saved at ${path}`))
