/**
 * The most minimal `task` implementation and its responsibilties.
 *
 * An implementation of `bionode-waterwheel` shall expose a user facing function
 * called `task` which takes two parameters:
 * - the first, `props`, is a plain object describing properties of the task
 * - the second, `operationCreator`, is a function which receives one parameter:
 *   a plain object called `opsProps` which describes the properties of the
 *   operation. `opsProps` **must not be referentially equal to `props`**.
 * - returns a function `invocableTask` which will accept a callback `cb` and a
 *   context `ctx`. the callback shall be called on completion or error of the
 *   `operation`.
 */

const _ = require('lodash')
const Promise = require('bluebird')

const { defaultCtx } = require('../lib/ctx')
const { generateUid } = require('./lifecycle')
const { createTask } = require('./reducers/tasks.js')

/**
 * The task factory function.
 *
 * @param props {Object} properties of this task
 * @param operationCreator {Function} function to be called with props
 */
const task = (store) => (props, operationCreator) => {
  // Type checking on params
  if (!_.isPlainObject(props) || !_.isFunction(operationCreator)) {
    throw new Error('incorrect parameter types: task(props, operationCreator) where props is a plain object and operationCreator is a function')
  }


  // Returns an invocable task that can be passed a cb and optionally a ctx
  // By defauult cb is a noop and ctx is defaultCtx
  // Provides promise and callback to user
  const invocableTask = (cb = _.noop, ctx = defaultCtx) => new Promise((resolve, reject) => {
    // First thing we need to do is create a `uid` for this task
    // The `uid` is a hash of the `input`, `output`, and `params` hashes
    // `hashes` is an object with input, output, params keys
    const { uid, hashes } = generateUid(props)

    // This kicks off the task lifecycle saga
    // operationCreator, taskResolve, and taskReject
    // will not be used in the reducer, BUT, they will be captured by the saga
    store.dispatch(createTask({
      uid,
      hashes,
      props,
      operationCreator,
      taskResolve: resolve,
      taskReject: reject
    }))

    let taskState

    // Case 1: resumable === 'on'
    //  Case 1A: output can be resolved
    //    1. create
    //    2. checkResumable
    //    3. afterOperation
    //    4. afterResolveOutput
    //  Case 1B: output cannot be resolved
    //  TODO fix infinite loop before 5/3/4
    //    1. create
    //    2. checkResumable
    //    3. afterOperation
    //    4. beforeOperation
    //    5. afterOperation
    //    6. afterResolveOutput
    // Case 2: resumable === 'off'
    //  1. create
    //  2. checkResumable
    //  3. beforeOperation
    //  4. afterOperation
    //  5. afterResolveOutput

    // Create the task
    // console.log('Creating the task')
    // lifecycle.create(props, ctx)
    //   .then((results) => {
    //     taskState = results

    //     console.log('Checking if resumable')
    //     return lifecycle.resumable(taskState)
    //   })
    //   // Check resumability: go to afterOperation or beforeOperation
    //   .then(isResumable => isResumable ? afterOperation(taskState) : beforeOperation(taskState))
    //   .catch(err => console.error(err))

    // Was either not resumable or resumable && resolveOutput failed
    function beforeOperation (taskState) {
      lifecycle.resolveInput(taskState)
      .then((results) => {
        Object.assign(taskState, { resolvedInput: results.resolvedInput })

        return lifecycle.createOperation(taskState, operationCreator)
      })
      .then((results) => {
        const { operation } = results
        taskState.dir = results.dir

        lifecycle.settleOperation(operation)
          .then((results) => {
            afterOperation(taskState)
          })
          .catch(err => console.error(err))
      })
      .catch(err => console.error(err))
    }

    // Either just after resumable check or after settleOperation
    function afterOperation(taskState) {
      console.log('Resolving output')
      lifecycle.resolveOutput(taskState)
        .then((results) => {
          Object.assign(taskState, results)

          afterResolveOutput(taskState)
        })
        .catch((err) => {
          console.error(err)
          beforeOperation(taskState)
        })
    }

    function afterResolveOutput (taskState) {
      lifecycle.validateOutput(taskState)
        .then((results) => {
          Object.assign(taskState, results)

          return lifecycle.postValidation(taskState)
        })
        .then((results) => {
          Object.assign(taskState, results)

          finish(taskState)
        })
        .catch(err => console.log(err))
    }

    function finish (taskState) {
      resolve(taskState)
    }
  }).asCallback(cb)

  return invocableTask
}

module.exports = task
