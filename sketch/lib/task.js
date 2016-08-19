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
 * Internally:
 * @param dispatch {Function} redux store.dispatch passed from above
 *
 * User facing:
 * @param props {Object} properties of this task
 * @param operationCreator {Function} function to be called with props
 * @returns {Function}
 */
const task = (dispatch) => (props, operationCreator) => {
  // Type checking on params
  if (!_.isPlainObject(props) || !_.isFunction(operationCreator)) {
    throw new Error('incorrect parameter types: task(props, operationCreator) where props is a plain object and operationCreator is a function')
  }

  // First thing we need to do is create a `uid` for this task
  // `hashes` is an object with input, output, params keys
  // The `uid` is a hash of the `input`, `output`, and `params` hashes
  const { uid, hashes } = generateUid(props)
  // console.log('params for uid before invoked', uid, props.params)

  // Returns an invocable task that can be passed a callback (or resolved as a
  // promise) and/or context
  const invocableTask = (cb = _.noop, ctx = defaultCtx) => new Promise((resolve, reject) => {
    // This kicks off the task lifecycle saga. Lifecycle saga will continue from
    // point just after this action has update the store.
    dispatch(createTask({
      uid,
      hashes,
      props,
      // For saga but not reducer
      operationCreator,
      taskResolve: resolve,
      taskReject: reject
    }))
  }).asCallback(cb)

  // Attach info to returned invocable task
  // Used to let join() know if it should stream across tasks
  // Could be used for a dry run using just glob patterns
  invocableTask.info = Object.assign(props, { uid })

  return invocableTask
}

module.exports = task
