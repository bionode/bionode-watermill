'use strict'

const globby = require('globby')

const types = require('./types.js')
const wrappers = require('./wrappers.js')

const { File } = types

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
  } else if (input instanceof Array && input[0] instanceof File) {
    input = input.map(file => file.value)
    globby(input).then((data) => {
      console.log(data)
      props.input = data
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

const join = function(...tasks) {
  // hack to work for parallel
  if (tasks[0] instanceof Array) {
    tasks = tasks[0]
  }

  // return tasks.reduceRight((prev, curr) => curr(prev()))
  let result = tasks[tasks.length-1]()
  for (let i=tasks.length-2; i>=0; i--) {
    result = tasks[i](result)
  }
  return result
}

const parallel = ({ taskLists, next }) => {
  let max = taskLists.length
  let count = 0

  function after() {
    return task(
    {
      input: null,
      output: null
    },
    () => {
      count++
      if (count === max) {
        console.log('PARALLEL FINISHED')
        next()()
      }
    }
    )()
  }

  const joinedTasks = []

  taskLists.forEach(taskList => {
    taskList.push(after)
    joinedTasks.push(join(taskList))
  })

  // Start parallel tasks
  return () => {
    joinedTasks.forEach(task => task())
  }
}

module.exports = {
  task,
  join,
  parallel,
  types,
  wrappers
}
