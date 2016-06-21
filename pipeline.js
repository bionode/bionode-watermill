'use strict'

// === WATERWHEEL ===
const { spawn } = require('child_process')
const globby = require('globby')

const task = (props, cb, next) => {
  const { output } = props

  const stream = cb(props)

  // Check current dir for files matching the pattern provided to the task
  // Then call next task in the pipeline with the resolved files
  stream.on('end', () => globby(output).then(next))

  return stream
}

const shell = (cmd) => {
  console.log('Starting: ' + cmd)
  cmd = cmd.split(' ')

  const process = spawn(cmd[0], cmd.slice(1))

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

// === PIPELINE ===
const ncbi = require('bionode-ncbi')

const config = {
  sraAccession: '2492428'
}

// thunk'it here or inside task?
// const/let does not get hoisted, and it is unnatural to describe pipelines
// in reverse order:
// const task2 = task(...)
// const task1 = task(,,task2)
// functions do get hoisted, but then we have a somewhat less pretty indented
// return task(...) boilerplate

// const samples = task(...)
// const samples = () => task(...)
function samples() {
  // task(props, cb, next)
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
    ( ({ input }) => ncbi.download(input.db, input.accession) ),
    // the next task for this task
    // data is the resolved output
    (data) => shell(`fastq-dump --split-files --skip-technical --gzip ${data[0]}`)
  )
}

// wip
// function dumps() {
//   return task({
//     input: '*.sra',
//     output: [1, 2].map(n => `*_${n}.fastq.gz`)
//   }, )
// }

const pipeline = samples

pipeline()
  .on('data', console.log)
  .on('end', () => console.log('Finished SRA download'))
