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
  // 1. Creating -> new entry in store w/ defaults + input/output/params
  // 2. is resumable -> **on** or **off** -> set status in store
  // 3. resolve input -> set task.resolvedInput
  //  - from props <- lineage
  //  - from fs
  // 4. action = actionCreator(resolvedInput) -> {process, promise, stream, curry cb} -> set task.action
  // 5. apply action to readable/writable of duplex
  // 6. resolve output -> set task.resovledOutput
  //  - traverse over output, over validators
  // 7. catch end/finish/close and destroy duplex, set resolvedOutput
  // Each task action returns a promise that resolves to the task uid.
  // ========================================================================
  // 1. Creating
  // ========================================================================
  actions.createTask(props)
  .then(uid => actions.checkResumable(uid))
  // ========================================================================
  // 2. Finding out if we need to validate for resumable task
  // ========================================================================
  .then(function (uid) {
    let { status } = getTask(uid)

    switch (status) {
      case CHECKING_RESUMABLE:
        // ========================================================================
        // 2a. If so, resolve output
        // ========================================================================
        return actions.resolveOutput(getTask(uid))
        break
      case RESOLVING_INPUT:
        // ========================================================================
        // 2b. Else, continue as normal
        // ========================================================================
        console.log('Skip to resolving input')
        throw new Error('TODO, resolving input')
        break
      default:
        throw new Error('status was not CHECKING_RESUMABLE or RESOLVING_INPUT')
    }
  })
  .then((uid) => {
    if (getTask(uid).resolvedOutput) {
      // ========================================================================
      // 2a -> FINISHED
      // ========================================================================
      stream._output = getTask(uid).resolvedOutput
      stream.destroy()

      throw new Error('EARLY STOP')
    } else {
      // ========================================================================
      // 2b -> 3 Resolve Input
      // ========================================================================
      return actions.resolveInput(getTask(uid))
    }
  })
  .then(function (uid) {
    // ========================================================================
    // 4. creating the action
    // ========================================================================
    return actions.createAction(getTask(uid), actionCreator)
  })
  .then(function (uid) {
    // ========================================================================
    // 5. set writable and/or readable of Duplex from action
    // ========================================================================

    // TODO add stream to store?
    return actions.setDuplex(getTask(uid), stream)
  })
  .then(function (uid) {
    // ========================================================================
    // 6. resolve output
    // ========================================================================
    return actions.resolveOutput2(getTask(uid), stream)
  })
  .then(function (uid) {
    // ========================================================================
    // 7. catch end/finish/close and destroy duplex, make output available
    // ========================================================================
    return actions.catchTask(stream)
  })
  .catch(err => console.log(err))

  return stream
}

module.exports = task
