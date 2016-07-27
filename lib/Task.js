'use strict'

const duplexify = require('duplexify')

const {
  CHECKING_RESUMABLE,
  RESOLVING_INPUT
} = require('./constants/task-status-types.js')

/**
 * The Task class.
 * @param {Object} store the redux store for this waterwheel instance
 * @params {Object} actions list of redux actions bound to dispatch passed in
 * internally from waterwheel. see reducers/task(s).js for these actions
 *
 * Public API:
 * @param  {Object} props         object with input and output
 * @param  {function} actionCreator function to produce an actionCreator
 * @return {stream}               stream to be orchestrated
 */
const task = (store, actions) => (props, actionCreator) => (passedOutput, trajectory) => {
  let stream = duplexify(null, null)

  const getTask = (uid) => store.getState().tasks[uid]

  // Run through the task lifecycle.
  // 1. Creating -> new entry in store w/ defaults + input/output/params
  // 2. is resumable -> **on** or **off** -> set status in store -> skip to 6
  // 3. resolve input -> set task.resolvedInput
  //  - from props <- lineage
  //  - from fs
  // 4. action = actionCreator(resolvedInput) -> {process, promise, stream, curry cb} -> set task.action
  // 5. set writable and/or readable of Duplex from action
  // 6. resolve output -> set task.resolvedOutput
  //  - traverse over output, over validators
  // 7. catch end/finish/close and destroy duplex, set resolvedOutput
  // Each task action returns a promise that resolves to the task uid.
  actions.createTask(props)
  .then(uid => actions.checkResumable(uid))
  .then((uid) => {
    const { status } = getTask(uid)

    switch (status) {
      case CHECKING_RESUMABLE:
        return actions.resolveOutput(getTask(uid))
        break
      case RESOLVING_INPUT:
        return Promise.resolve(uid)
        // TODO resolving input
        break
      default:
        throw new Error('status was not CHECKING_RESUMABLE or RESOLVING_INPUT')
    }
  })
  .then((uid) => {
    // TODO not this hack
    stream.stateWormhole = getTask(uid)

    return Promise.resolve(uid)
  })
  .then((uid) => getTask(uid).resolvedOutput ? actions.runValidators(getTask(uid), 'before') : Promise.resolve(uid))
  .then((uid) => {
    if (getTask(uid).validated) {
      finish(uid)
    } else {
      actions.resolveInput(getTask(uid), passedOutput, trajectory)
      .then(uid => actions.createAction(getTask(uid), actionCreator))
      .then(uid => actions.setDuplex(getTask(uid), stream))
      .then(uid => actions.catchTask(getTask(uid), stream))
      .then(uid => actions.resolveOutput(getTask(uid)))
      .then(uid => actions.runValidators(getTask(uid), 'after'))
      .then((uid) => {
        if (getTask(uid).validated) {
          finish(uid)
        }
      })
      .catch(err => console.log(err))
    }
  })
  .catch(err => console.log(err))

  function finish(uid) {
    stream._output = getTask(uid).resolvedOutput
    stream.destroy()
    // TODO not this
    stream.emit('destroy')
  }

  return stream
}

module.exports = task
