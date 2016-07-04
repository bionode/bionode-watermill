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

const traverser = require('./traverser.js')
const { tab } = require('./utils.js')
// TODO refactor this out
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
  let resolvedInput
  if (passedOutput) {
    // a. use values passed in from last Task OR
    resolvedInput = traverser(passedOutput, (key, obj) => {
      const val = obj[key]

      switch(key) {
        case 'value':
          obj = val
          console.log(tab(2) + 'resolved ' + chalk.magenta(val) + ' from ' + chalk.cyan('value'))
          break
        case 'file':
          // Compare to task input
          // Can also check obj.resolved here
          // TODO handle better input object 
          const matches = multimatch([val], [input.file])

          if (matches.length === 1) {
            obj = val
            console.log(tab(2) + 'input: ' + chalk.magenta(val) + ' from passed in ' + chalk.cyan('file'))
          }
          break
        case 'resolved':
          break
        default:
          console.log(key, val)
          return false
      }

      return obj
    })
  } else {
    // b. resolve props to existing files
    resolvedInput = traverser(input, (key, obj) => {
      const val = obj[key]

      switch(key) {
        case 'value':
          obj = val
          console.log(tab(2) + 'resolved ' + chalk.magenta(val) + ' from ' + chalk.cyan('value'))
          break
        case 'file':
          // Compare to task input
          // Can also check obj.resolved here
          // TODO handle better input object 
          const matches = globby.sync(val)

          if (matches.length >= 0) {
            obj = matches[0]
            console.log(tab(2) + 'input: ' + chalk.magenta(matches[0]) + ' from ' + chalk.cyan('file'))
          }
      }


      return obj
    })
  }
  
  

  // 3. creating
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

    return traverser(output, (key, obj) => { 
      const val = obj[key]

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
          return false
      }

      return obj
    })
   
    // TODO check for resolved:true
  }

  // Function to be called from outside
  stream.output = resolveOutput

  if (isWritable(stream) && !isReadable(stream)) {
    stream.on('finish', function (chunk)  {
      console.log(tab(1) + 'Stream has finished')

      if (chunk) {
        console.log('Last chunk ' + chunk)
      } 
     
      // this.output = resolveOutput()
    }) 
  } else if (stream.stdin && stream.stdout && stream.stderr) {
    stream.on('close', () => {
      console.log(tab(1) + 'Child process finished')
      
      // stream.output = resolveOutput()
    })
  }
  // 6. ending, passing on output
  
  // console.log(stream)

  return stream
}

module.exports = task
