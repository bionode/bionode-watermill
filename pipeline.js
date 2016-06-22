'use strict'

// === WATERWHEEL ===
const { spawn } = require('child_process')
const globby = require('globby')

// mostly for instanceof checks inside task
class File {
  constructor(value) {
    this.value = value
  }
}

const task = (props, cb) => (next) => () => {
  console.log('Starting a task')
  let { input, output } = props

  // TODO recursively traverse object
  // TODO use generators to yield Promises to be cleaner
  if (input instanceof File) {
    // If we got special classes as input, they are either stdin or files,
    // so we need to resolve the pattern
    globby(input.value)
      .then((data) => {
        // TODO unhardcoded to one file
        props.input = data[0]
        return ready()
      })
  } else {
    return ready()
  }

  function ready() {
    const stream = cb(props)

    let ended = false

    const handleFinish = () => {
      if (!ended) {
        ended = true
      } else {
        // Don't start shell process more than once
        return null
      }

      let length = 1

      if (typeof(output) === 'string') {
        length = 1
      } else if (output instanceof File) {
        output = output.value
        length = 1
      } else if (output instanceof Array && output[0] instanceof File) {
        output = output.map(file => file.value)
        length = output.length
      }

      globby(output).then((data) => {
        console.log('Wanted: ')
        console.log(output)
        console.log(data.length === length ? 'File(s) produced as expected: ' : 'File(s) not created: ')
        console.log(data)
        next()
      })
    }

    // Check current dir for files matching the pattern provided to the task
    // Then call next task in the pipeline with the resolved files
    stream.on('end', handleFinish)
    // fs.createWriteStream throws finish
    stream.on('finish', handleFinish)
    // bash processes
    stream.on('close', handleFinish)

    return stream
  }
}

function join(...tasks) {
  // return tasks.reduceRight((prev, curr) => curr(prev()))
  let result = tasks[tasks.length-1]()
  for (let i=tasks.length-2; i>=0; i--) {
    result = tasks[i](result)
  }
  return result

  // return tasks[0](tasks[1](tasks[2]()))
  // // return tasks.reduce( (prev, curr) => prev(curr()) )
  // return tasks.reduce((prev, curr, index, arr) => {
  //   console.log(prev)
  //   console.log(curr)
  //   console.log(index)
  //   return prev(curr())
  // })
  // combined({}, () => console.log('after'))
}

const shell = (cmd) => {
  console.log('Starting: ' + cmd)
  cmd = cmd.split(' ')

  const process = spawn(cmd[0], cmd.slice(1))

  process.stdout.on('data', (data) => console.log(`stdout: ${data}`) )

  process.stderr.on('data', (data) => console.log(`stderr: ${data}`) )

  process.on('close', (code) => console.log(`child process exited with code ${code}`) )

  return process
}

// === UTILS ===
const last = (str, sep) => {
  const splitted = str.split(sep)
  return splitted[splitted.length - 1]
}

// === PIPELINE ===
const fs = require('fs')
const ncbi = require('bionode-ncbi')
const request = require('request')

const config = {
  sraAccession: '2492428',
  referenceURL: 'http://ftp.ncbi.nlm.nih.gov/genomes/all/GCA_000988525.2_ASM98852v2/GCA_000988525.2_ASM98852v2_genomic.fna.gz'
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

function after(next) {
  return task(
  {
    input: null,
    output: null
  },
  () => shell('echo AFTER')
  )(next)
}


// join will populate the nexts
const pipeline = join(samples, fastqDump, after)

pipeline()
//   .on('data', console.log)
//   .on('end', () => console.log('Finished SRA download'))


// const downloadAndIndex = join(downloadReference, bwaIndex, after)
// downloadAndIndex()
  // .on('close', () => console.log('sdsds'))
