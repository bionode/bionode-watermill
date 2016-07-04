'use strict'

// Can potentially be replaced by explicty type objects
const isStream = require('isstream')
const { isWritable, isReadable, isDuplex } = isStream
// stream-from-x modules
const StreamFromPromise = require('stream-from-promise')
// TODO make a module of this
const StreamFromCallback = require('./stream-from-callback.js')
// pattern matching
const globby = require('globby')
const multimatch = require('multimatch')

const chalk = require('chalk')

const resolveObj = require('./resolveObj.js')
const { tab } = require('./utils.js')
// TODO refactor this out
const { File } = require('./types.js')
const { Process } = require('./wrappers.js')


/**
 * The Task class.
 * @param  {Object} props         object with input and output
 * @param  {function} actionCreator function to produce an actionCreator
 * @return {stream}               stream to be orchestrated
 */
const task = (props, actionCreator) => (passedOutput) => {
  // TODO check if actionCreator is a function
  if (props === null) {
    console.log('Got null props')
    props = {
      input: null,
      name: 'noname'
    }
  } 


  // TODO handle unnamed task
  const { input, output, name } = props 
  // TODO task lifecycle
  // 1. preparing
  console.log(`Running task ${chalk.blue(name)}`)
  // 2. resolving input
  console.log(tab(1) + chalk.italic('Resolving input'))
  
  // TODO handle passed in input
  // if (passedOutput) {
  //   // a. use values passed in from last Task OR
  //   objectWalk(passedOutput)
  // }
  
  const resolvedInput = resolveObj (input)

  // console.log('Before:')
  // console.log(input)
  // console.log('After: ')
  // console.log(resolvedInput)

  // if (passedOutput) {
  //   // a. use values passed in from last Task OR
  //   if (passedOutput instanceof File && input instanceof File) {
  //     const matches = multimatch([passedOutput.value], [input.value])
  //     
  //     if (matches.length === 1) {
  //       resolvedInput = { input: matches[0] }
  //       console.log(tab(2) + 'input: ' + chalk.magenta(resolvedInput.input) + ' from passed in ' + chalk.cyan('File'))
  //     }
  //   } 
  // } else {
  //   // b. resolve props to existing files
  //   if (input instanceof File) {
  //     try {
  //       const globbyData = globby.sync(input.value)
  //       resolvedInput = { input: globbyData[0] } 
  //     } catch(e) {
  //       console.error(e)
  //     }
  //     console.log(tab(2) + 'input: ' + resolvedInput.input + ' from ' + input.value)
  //   } else if (typeof(input) === 'string') {
  //     console.log(tab(2) + 'input: ' + chalk.magenta(resolvedInput.input) + ' from ' + chalk.cyan('value'))
  //   }
  //
  //   if (typeof(input) === 'object') {
  //     console.log(tab(2) + 'input:')
  //     for (let item in input) {
  //       if (typeof(input[item]) === 'string') {
  //         console.log(tab(3) + item + ': ' + input[item] + ' from ' + chalk.cyan('value'))
  //       }
  //     }
  //   }
  // }
  // 3. creating
  // TODO use input/output, resolved promise, callback to determine stream type.
  // Or just always be duplex.
  let stream
  // const stream = new Readable()

  // Create the action
  // TODO handle actual callbacks, in which case we cant call this yet 
  const resolvedProps = Object.assign(props, { input: resolvedInput }) 
  const action = actionCreator(resolvedProps)


  // Get its type: can be stream, promise, function
  const actionType = (() => {
    if (action instanceof Process) {
      return 'process'
    } else if (action instanceof Promise) {
      return 'promise'
    } else if (typeof(action) === 'function') {
      // assume action returns function of signature fn(err, data)
      return 'curried-callback'
    } else if (isStream(action)) {
      return 'stream'
    } else {
      return 'unkown'
    }
  })()

  // 4. running
  switch(actionType) {
    case 'process':
      stream = action.getSpawned()
      console.log(tab(1) + 'Creating a child process') 
      break
    case 'promise':
      stream = new StreamFromPromise(action)
      console.log(tab(1) + 'Creating a Readable from promise')
      break
    case 'curried-callback':
      stream = new StreamFromCallback(action)
      break
    case 'stream':
      stream = action
      console.log(tab(1) + 'Creating a stream from stream')
      break
    default:
      // TODO this is a hack to make join work for now
      stream = new Duplex()
      console.log('Bad action type')
  }

  // 5. resolving output
    
  const resolveOutput = () => {
    const resolveFile = ({ file }) => {
      console.log(tab(2) + chalk.yellow('Looking') + ' for ' + file)
      const globbyData = globby.sync(file)
      
      if (globbyData.length > 0) {
        console.log(tab(3) + chalk.green(' found ') + globbyData[0])
        return { file: globbyData[0], resolved: true }
      }

      return Object.assign(file, { resolved: false })
    }
  

    console.log(tab(1) + chalk.italic('Resolving output'))

    const traverse = (obj) => {
      if (obj instanceof Array) {
        return obj.map(item => resolveObj(item))
      }

      for (let key in obj) {
        let val = obj[key]

        switch (key) {
          case 'value':
            obj = val
            break
          case 'file':
            obj = resolveFile(obj)
            break
          // case 'fileResolved':
          //   console.log('got a fileResolved')
          //   break
          default:
            obj[key] = traverse(val)
        }
      }

      return obj
    }
    const resolvedOutput = traverse(output)

    // if (output instanceof Array) {
    //   if (output[0] instanceof File) {
    //     const resolvedFiles = output.map(file => resolveFile(file))
    //     
    //     const passes = resolvedFiles.every(file => file)
    //
    //     if (passes) {
    //       stream.emit('task.done', { file: output[0] })
    //     } else {
    //       // TODO this
    //       stream.emit('task.failedOutputResolution', output)
    //     }
    //   } 
    // } else if (output instanceof File) {
    //   resolveFile(output)
    //
    //   stream.emit('task.done', output)
    // }
  }

  if (isWritable(stream) && !isReadable(stream)) {
    stream.on('finish', (chunk) => {
      console.log(tab(1) + 'Stream has finished')
      if (chunk) {
        console.log('Last chunk ' + chunk)
      } 
     
      resolveOutput() 
    }) 
  } else if (stream.stdin && stream.stdout && stream.stderr) {
    stream.on('close', () => {
      console.log(tab(1) + 'Child process finished')

      resolveOutput()
    })
  }
  // 6. ending, passing on output
  
  // console.log(stream)

  return stream
}

module.exports = task
