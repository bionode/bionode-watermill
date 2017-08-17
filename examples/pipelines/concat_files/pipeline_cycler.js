'use strict'

// === WATERMILL ===
const {
  task,
  join,
  junction,
  fork
} = require('../../..')

const fs = require('fs')
const path = require('path')

// wrap every task and pipeline of a file in a single function (in this case
// 'pipeline' function
const fileFlow = (infile) => {
  const concatTask = task({
    // notice that the input file (infile) is here passed as a param of the
    // pipeline function
      input: infile,
      output: '*.fas',
      params: {output: 'concat.fas'}
    }, ( object ) => {
      console.log('input_directory: ', object.dir)
      console.log('input_variable: ', object, object.params)
      return `cat ${object.input} > ${object.params.output}`
    }
  )

  const task0 = task({name: 'coco'}, () => `echo "something0"`)

  const task1 = task({name: 'nut'}, () => `echo "something1"`)

  const task2 = task({name: 'foo'}, () => `echo "something2"`)

  // this returns the pipeline for a file
  const pipelineAssemblage = join(concatTask, task0, fork(task1, task2))
  return pipelineAssemblage
}

// checks for files in current working directory (where the script is executed)
fs.readdir(process.cwd(), (err, files) => {
  files.forEach(infile => {
    // checks file extension, in this case '.fas'
    if (path.extname(infile) === '.fas') {
      // executes the pipeline function for each infile
      const pipelineMaster = fileFlow(infile)
      pipelineMaster()
      // this is nice because in fact each pipelineMaster creates its own DAG
      // however graphson output file may be corrupt because DAGs are being
      // fired for each file, so last one to get set is the one saved to
      // graphson.json
      // TODO fix this graphson output file behavior
    }
  })
})