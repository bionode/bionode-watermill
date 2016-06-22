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
  let { input, output } = props

  // TODO recursively traverse object
  // TODO use generators to yield Promises to be cleaner
  if (input instanceof File) {
    // If we got special classes as input, they are either stdin or files,
    // so we need to resolve the pattern
    globby(input.value)
      .then((data) => {
        props.input = data[0]
        return ready()
      })
  } else {
    return ready()
  }

  function ready() {
    const stream = cb(props)

    // Check current dir for files matching the pattern provided to the task
    // Then call next task in the pipeline with the resolved files
    stream.on('end', () => globby(output).then(next))

    return stream
  }
}

function join(...tasks) {
  // return tasks[0](tasks[1]())
  return tasks.reduce( (prev, curr) => prev(curr()) )
}

const shell = (cmd) => {
  console.log('Starting: ' + cmd)
  cmd = cmd.split(' ')

  const process = spawn(cmd[0], cmd.slice(1))

  process.stdout.on('data', (data) => console.log(`stdout: ${data}`) )

  process.stderr.on('data', (data) => console.log(`stderr: ${data}`) )

  process.on('close', (code) => console.log(`child process exited with code ${code}`) )
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

function fastqDump() {
  return task({
    input: new File('**/*.sra'),
    output: [1, 2].map(n => `*_${n}.fastq.gz`)
  },
  ({ input }) => shell(`fastq-dump --split-files --skip-technical --gzip ${input}`)
  )()
}

// pass in next which is given globbied output
// const pipeline = samples( (data) => shell(`fastq-dump --split-files --skip-technical --gzip ${data[0]}`) )
// or nothing to stop the pipeline
// const pipeline = samples()

// join will populate the nexts
const pipeline = join(samples, fastqDump)

pipeline()
  .on('data', console.log)
  .on('end', () => console.log('Finished SRA download'))
