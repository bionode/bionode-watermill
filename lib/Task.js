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

const {
  CHECKING_RESUMABLE,
  RESOLVING_INPUT
} = require('./constants/task-status-types.js')

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
const task = (store, actions) => (props, actionCreator) => (passedOutput) => {
  let stream = duplexify(null, null)

  const getTask = (uid) => store.getState().tasks[uid]

  // Run through the task lifecycle. 
  // Each task action returns a promise that resolves to the task uid.
  // ========================================================================
  // 1. Creating
  // ========================================================================
  actions.createTask(props)
  .then(function (uid) {
    console.log('redux: ', `Running task ${chalk.blue(getTask(uid).name)}`)

    setTimeout(() => stream.emit('wow', uid), 1000)

    return actions.checkResumable(uid)
    // return new Promise((resolve, reject) => resolve())
  })
  // ========================================================================
  // 2. Finding out if we need to validate for resumable task
  // ========================================================================
  .then(function (uid) {
    let { status } = getTask(uid)

    switch (status) {
      case CHECKING_RESUMABLE:
        // ========================================================================
        // 3a. If so, resolve output
        // ========================================================================
        console.log(tab(1) + 'redux: Checking if valid output already exists')
        return actions.resolveOutput(getTask(uid))
        break
      case RESOLVING_INPUT:
        // ========================================================================
        // 3b. Else, continue as normal
        // ========================================================================
        console.log('Skip to resolving input')
        break
      default:
        throw new Error('status was not CHECKING_RESUMABLE or RESOLVING_INPUT')
    }
  })
  .then((uid) => {
    if (getTask(uid).resolvedOutput) {
      console.log('redux: Could resolve output before running')
      // stream.output = () => console.log('resolvedOutput: ', getTask(uid).resolvedOutput)
      stream._output = getTask(uid).resolvedOutput
      stream.destroy()
    } else {
      console.log('redux: Could not resolve output before running')

      // stubCloser()
      newRest()

      // Something happens here that stops .on('close') being caught from
      // outside
      // rest()

      // stream.on('close', function() {
      //   console.log('Internally received close')
      //   stream.destroy()
      // })

      // TODO next, resolveInput
    }
  })
  .catch(err => console.log(err))

  function stubCloser () {
    setTimeout(() => {
      console.log('Deferred close')
      stream.output = () => console.log('STUB OUTPUT')
      stream.destroy()
    }, 1000)
  }

  function newRest() {
    const { input, output, name } = props
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
      console.log('resolving to fs')
      resolvedInput = outerResolver(input, props, 'start')
    }

    const resolvedProps = Object.assign(props, { input: resolvedInput })
    console.log('resolvedProps', resolvedProps)
    const action = actionCreator(resolvedProps)

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

    console.log(actionType)

    if (isReadable(action) && !isWritable(action)) {
      stream.setReadable(action)
      console.log('isReadable')
    } else if (isWritable(action) && !isReadable(action)) {
      console.log('isWritable')
      stream.setWritable(action)
    }

    stream.on('finish', () => console.log('internally received finish'))

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

    // IDEA: create stream._output on end, finish, events, then make
    // a close by calling stream.destroy()
    // (but wont work for child process)?

    stream.on('finish', function (chunk)  {
      console.log(tab(1) + 'Stream has finished')

      stream._output = stream.output()
      stream.destroy()

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

    // setTimeout(() => {
    //   console.log('Deferred close')
    //   stream.output = () => console.log('STUB OUTPUT')
    //   stream.destroy()
    // }, 1000)
  }

  // TODO check if actionCreator is a function

  // rest()
  function rest () {
    const { input, output, name } = props
    // 2. Checking if task skippable
    if (false) {
      console.log(tab(1) + 'first: Checking if valid output already exists')
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

        console.log('GONNA DESTROY STREAM')
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
        console.log('elsed to make a duplex')
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
  }

  // 7. ending, passing on output
  return stream
}

module.exports = task
