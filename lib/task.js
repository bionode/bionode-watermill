/**
 * The most minimal `task` implementation and its responsibilties.
 *
 * An implementation of `bionode-watermill` shall expose a user facing function
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

const { defaultContext } = require('../lib/constants/default-task-state.js')
const { generateUid } = require('./lifecycle')
const { createTask } = require('./reducers/tasks.js')

const hash = require('./utils/hash.js')
const stringify = require('json-stable-stringify')

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
const task = (dispatch) => (props, operationCreator, instanceUid = null) => {
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
  // const invocableTask = (cb = _.noop, ctx = defaultContext) => {
  const invocableTask = (cb = _.noop, ctx = defaultContext) => {
    // const self = this
    return new Promise((resolve, reject) => {
      // Take trajectory from context
      const { trajectory } = ctx
      // Make a new context for the task we are about to create using it
      const context = { trajectory }

      // creates new uid with merkle tree like hash
      // Although task can only have one uid, here we assure that it can be
      // duplicated in a pipeline since it needs the previous context to
      // generate newUid
      const newUid = childHashesToUid(context, uid)

      // This kicks off the task lifecycle saga. Lifecycle saga will continue from
      // point just after this action has update the store.
      dispatch(createTask({
        uid: instanceUid ? instanceUid : newUid,
        hashes,
        props,
        context,
        // For saga but not reducer
        operationCreator,
        taskResolve: resolve,
        taskReject: reject
      }))
    }).asCallback(cb)
  }
  // const invocableTask = (cb = _.noop, ctx = defaultContext) => new Promise((resolve, reject) => {

  // Attach info to returned invocable task
  // Used to let join() know if it should stream across tasks
  // Could be used for a dry run using just glob patterns

  // const obj = {
  //   instanceUid: null
  // }
  // invocableTask = invocableTask.bind(obj)

  invocableTask.info = Object.assign({
    uid: instanceUid ? instanceUid : uid,
    props,
    operationCreator,
    type: 'task'
  })
  // invocableTask.setInstanceUid = (instanceUid) => {
  //   obj.instanceUid = instanceUid
  //   invocableTask.info.uid = instanceUid
  // }

  return invocableTask
}

// function to generate new hash with for merkle tree
const childHashesToUid = (context, uid) => {
  // TODO some instances seem to return undefined (not sure which ones)
  if (context.trajectory !== undefined) {
    const childHashes = Object.keys(context).reduce(
      (sum, key) => Object.assign({}, sum, {[key]: hash(stringify(context[key]))}),
      {}
    )
    uid = hash(childHashes.trajectory + uid)
    return uid
  } //else {
    //console.log('undefined trajectory: ', context, uid)
  //}
}

module.exports = task
