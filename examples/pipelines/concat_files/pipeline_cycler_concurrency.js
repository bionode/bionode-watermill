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
const Promise = require('bluebird')

// makes concurrency the third argument of code in shell
const concurrency = parseFloat(process.argv[2] || "Infinity")

/**************************/
/* PIPELINE FOR EACH FILE */
/**************************/

// wrap every task and pipeline of a file in a single function (in this case
// 'pipeline' function
const fileFlow = (infile) => {
  console.log('inputfile', infile)
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

/************************************************/
/* CONTROL THE FLOW OF INPUTS INTO THE PIPELINE */
/************************************************/

// checks for files in current working directory (where the script is executed)
fs.readdir(process.cwd(), (err, files) => {
  const desiredFiles = []
  files.forEach(infile => {
    // checks file extension, in this case '.fas' and creates an array with
    // just these files
    if (path.extname(infile) === '.fas') {
      desiredFiles.push(infile)
    }
  })
  // bluebird Promise.map is used here for concurrency where we can
  // limit the flow of files into the pipeline. So, if concurrency is set
  // to 1 it will first execute the first file, then the second and so on.
  // But if 2 it will execute the first 2 files and then the other 2.
  return Promise.map(desiredFiles, (infile) => {
    const pipelineMaster = fileFlow(infile)
    return pipelineMaster()
  }, {concurrency})
})