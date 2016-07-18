'use strict'

// Can potentially be replaced by explicty type objects
const isStream = require('isstream')
const { isWritable, isReadable, isDuplex } = isStream
const duplexify = require('duplexify')
// stream-from-x modules
const StreamFromPromise = require('stream-from-promise')
// TODO make a module of this
const StreamFromCallback = require('./stream-from-callback.js')


const chalk = require('chalk')

const { tab, isChildProcess } = require('./utils.js')
// TODO refactor this out
const { Process } = require('./wrappers.js')
const outerResolver = require('./resolvers.js')
const validators = require('./validators.js')

/**
 * The Task class.
 * @params {Object} actions list of redux actions bound to dispatch passed in
 * internally from waterwheel.
 *
 * Public API:
 * @param  {Object} props         object with input and output
 * @param  {function} actionCreator function to produce an actionCreator
 * @return {stream}               stream to be orchestrated
 */
const task = (actions) => (props, actionCreator) => (passedOutput) => {

  let stream = duplexify(null, null)

  actions.addTask(props).then(uid => {                          
    console.log('got uid: ', uid)
    
    
    // TODO check if actionCreator is a function
    if (props === null) {
      console.log('Got null props')
      props = {
        input: null,
        name: 'noname'
      }
    }
    
    // Task props defaults
    props = Object.assign({}, { skippable: true }, props)


    // TODO handle unnamed task
    const { input, output, name } = props 
    // TODO task lifecycle
    // 1. preparing
    console.log(`Running task ${chalk.blue(name)}`)
    // 2. Checking if task skippable
    if (output && !output.stream && props.skippable) {
      console.log(tab(1) + 'Checking if valid output already exists')
      let alreadyOutput = outerResolver(output, props, 'checking', { 
        fileResolvers: [validators.matchHash, validators.matchTime]
      })
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
        console.log(tab(1) + 'Skipping this task')

        stream.output = () => {
          return alreadyOutput
        }

        stream.destroy()

        return stream
      }
    }

    // 3. resolving input
    console.log(tab(1) + chalk.italic('Resolving input'))
    
    let resolvedInput
    if (passedOutput && !props.skipPassed) {
      // a. use values passed in from last Task OR
      resolvedInput = outerResolver(passedOutput, props, 'passed')
      console.log('Resolved from state:')
      console.log(resolvedInput)

      // pop all DELETE_ME's
      // may break for more than one DELETE_ME..
      resolvedInput.forEach((item, i) => item === 'DELETE_ME' ? resolvedInput.splice(i, 1) : null)

      console.log(resolvedInput)
    } else {
      // b. resolve props to existing files
      resolvedInput = outerResolver(input, props, 'start') 
    }
    

    // 4. creating
    // Create the action
    // TODO handle actual callbacks, in which case we cant call this yet 
    const resolvedProps = Object.assign(props, { input: resolvedInput }) 
    const action = actionCreator(resolvedProps)


    // Get its type: can be stream, promise, function
    const actionType = (() => {
      if (isChildProcess(action)) {
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
        stream.setReadable(action)
        console.log(tab(1) + 'Creating a child process') 
        break
      case 'promise':
        if (output && output.stream) {
          if (output.stream === 'object') {
            console.log(tab(1) + 'Creating a Readable for object from promise')

            // TODO do this earlier
            stream = duplexify.obj(null, null)
            stream.end()

            stream.setReadable(StreamFromPromise.obj(action))
          }
        } else {
          console.log(tab(1) + 'Creating a Readable for buffer/string from promise')
          stream.setReadable(StreamFromPromise(action))
        }
        break
      case 'curried-callback':
        stream.setReadable(new StreamFromCallback(action))
        break
      case 'stream':
        if (isReadable(stream) && !isWritable(stream)) {
          stream.setReadable(action)
        } else if (isWritable(stream) && !isReadable(stream)) {
          stream.setWritable(stream)
        } else {
          // Duplex
          stream = action
        }
        console.log(tab(1) + 'Creating a stream from stream')
        break
      default:
        // TODO this is a hack to make join work for now
        // stream = new Duplex()
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
     
      const customValidators = {
        fileResolvers: []
      }
      if (props.validators && props.validators.file) {
        customValidators.fileResolvers = props.validators.file
      }

      return outerResolver(output, props, 'end', customValidators)
    }

    // TODO handle this inside or outside task?
    // Wrap task with something that tells it if it is alone or not?
    stream.on('finish', function (chunk)  {
      console.log(tab(1) + 'Stream has finished')

      if (chunk) {
        console.log('Last chunk ' + chunk)
      }
    })

    stream.on('close', () => {
      console.log(tab(1) + 'Child process finished')
    })

    stream.on('data', function(data) {
      if (actionType === 'promise') {
        this._resolved = data
      }
    })
   
    stream.on('error', (err) => {
      console.log(err)
    })
  })

  // 7. ending, passing on output
  return stream
}

module.exports = task
