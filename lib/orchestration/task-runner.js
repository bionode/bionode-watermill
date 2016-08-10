'use strict'

const duplexify = require('duplexify')
const { EventEmitter2 } = require('eventemitter2')

const {
  CHECKING_RESUMABLE,
  RESOLVING_INPUT
} = require('../constants/task-status-types.js')

/**
 * The Task class.
 * @param {Object} store the redux store for this waterwheel instance
 * @params {Object} actions list of redux actions bound to dispatch passed in
 * internally from waterwheel. see reducers/task(s).js for these actions
 *
 * Public API:
 * @param  {Object} props         object with input and output
 * @param  {function} operationCreator function to produce an operationCreator
 * @return {stream}               stream to be orchestrated
 */
const task = (store, actions) => (props, operationCreator) => (trajection = []) => {
  // Use this to emit task.finish, task lifecycle events
  const emitter = new EventEmitter2()

  // Pull task from redux store
  const getTask = (uid) => store.getState().tasks[uid]

  // Run through the task lifecycle.
  // 1. Creating -> new entry in store w/ defaults + input/output/params
  // 2. is resumable -> **on** or **off** -> set status in store -> skip to 6
  // 3. resolve input -> set task.resolvedInput
  //  - from props <- lineage
  //  - from fs
  // 4. operation = operationCreator(resolvedInput) -> {process, promise, stream} -> set task.operation
  // 5. set writable and/or readable of Duplex from operation
  // 6. resolve output -> set task.resolvedOutput
  //  - traverse over output, over validators
  // 7. catch end/finish/close and destroy duplex, set resolvedOutput
  // Each task action returns a promise that resolves to the task uid.

  actions.createTask(props) // tasks/create
  .then(uid => actions.addTrajection(uid, trajection))
  .then(uid => actions.checkResumable(uid))
  .then((uid) => {
    // TODO make this check inside an action creator
    const { status } = getTask(uid)

    switch (status) {
      case CHECKING_RESUMABLE:
        return actions.resolveOutput(uid)
        break
      case RESOLVING_INPUT:
        return Promise.resolve(uid)
        // TODO resolving input
        break
      default:
        throw new Error('status was not CHECKING_RESUMABLE or RESOLVING_INPUT')
    }
  })
  .then((uid) => getTask(uid).resolvedOutput ? actions.runValidators(uid, 'before') : Promise.resolve(uid))
  .then((uid) => {
    if (getTask(uid).validated) {
      return Promise.resolve(uid)
    } else {
      return actions.resolveInput(uid, trajection)
      .then(uid => actions.createAction(uid, operationCreator))
      .then(uid => actions.resolveOutput(uid))
      .then(uid => actions.runValidators(uid, 'after'))
      .then((uid) => {
        if (getTask(uid).validated) {
          return Promise.resolve(uid)
        }
      })
      .catch(err => console.log(err))
    }
  })
  .then(uid => actions.addOutput(uid))
  .then(uid => finish(uid))
  .catch(err => console.log(err))

  function finish(uid) {
    emitter.emit('task.finish', getTask(uid))
  }

  emitter.async = new Promise((resolve, reject) => {
    emitter.on('task.finish', results => resolve(results))
  })

  return emitter
}

module.exports = task
