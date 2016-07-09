'use strict'

// Can potentially be replaced by explicty type objects
const isStream = require('isstream')
const { isWritable, isReadable, isDuplex } = isStream
// stream-from-x modules
const StreamFromPromise = require('stream-from-promise')
// TODO make a module of this
const StreamFromCallback = require('./stream-from-callback.js')


const chalk = require('chalk')

const { tab } = require('./utils.js')
// TODO refactor this out
const { Process } = require('./wrappers.js')
const outerResolver = require('./resolvers.js')

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
  // 2. Checking if task skippable
  if (output && !output.stream && props.skippable) {
    console.log(tab(1) + 'Checking if valid output already exists')
    let alreadyOutput = outerResolver(output, props, 'end')
    // Check that everything resolved
    let safe = true
    if (!(alreadyOutput instanceof Array)) {
      alreadyOutput = [alreadyOutput]
    }

    for (let item of alreadyOutput) {
      if (!item.resolved) {
        safe = false
      }
    }
    if (safe) {
      console.log('Can skip this task')
      let stream = require('readable-stream').Readable({
        read: function() {
          this.push(null)
        }
      })

      stream.output = () => {
        return alreadyOutput
      }

      stream.push(null)

      stream.on('data', function() {} )

      return stream
    }
  }

  // 3. resolving input
  console.log(tab(1) + chalk.italic('Resolving input'))
  
  let resolvedInput
  if (passedOutput) {
    // a. use values passed in from last Task OR
    resolvedInput = outerResolver(passedOutput, props, 'passed')
  } else {
    // b. resolve props to existing files
    resolvedInput = outerResolver(input, props, 'start') 
  }
  

  // 4. creating
  let stream
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

  // 5. running
  switch(actionType) {
    case 'process':
      stream = action.getSpawned()
      console.log(tab(1) + 'Creating a child process') 
      break
    case 'promise':
      if (output && output.stream) {
        if (output.stream === 'object') {
          console.log(tab(1) + 'Creating a Readable for object from promise')
          stream = new StreamFromPromise.obj(action)
        }
      } else {
        console.log(tab(1) + 'Creating a Readable for buffer/string from promise')
        stream = new StreamFromPromise(action)
        // stream.on('data', (data) => console.log('got data?: ', data.toString()))
      }
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

  // 6. resolving output
  
  // Called on completion of stream from join
  // TODO do this also for one task at a time 
  stream.output = () => {
    console.log(tab(1) + chalk.italic('Resolving output'))

    if (actionType === 'promise') {
      return stream._resolved
    }
    
    return outerResolver(output, props, 'end')
  }

  // TODO handle this inside or outside task?
  // Wrap task with something that tells it if it is alone or not?
  if (isWritable(stream) && !isReadable(stream)) {
    stream.on('finish', function (chunk)  {
      console.log(tab(1) + 'Stream has finished')

      if (chunk) {
        console.log('Last chunk ' + chunk)
      } 
    }) 
  } else if (stream.stdin && stream.stdout && stream.stderr) {
    stream.on('close', () => {
      console.log(tab(1) + 'Child process finished')
    })
  } else if (actionType === 'promise') {
    stream.on('data', function(data) {
      stream._resolved = data
    })
  }
 
  stream.on('error', (err) => {
    console.log(err)
  })

  // 7. ending, passing on output
  return stream
}

module.exports = task
