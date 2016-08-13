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
const asyncDone = require('async-done')

const lifecycleMethods = require('./lifecycle.js')

/**
 * The task factory function.
 *
 * @param props {Object} properties of this task
 * @param operationCreator {Function} function to be called with props
 */
const task = (props, operationCreator) => {
  // Type checking on params
  if (!_.isPlainObject(props) || !_.isFunction(operationCreator)) {
    throw new Error('incorrect parameter types: task(props, operationCreator) where props is a plain object and operationCreator is a function')
  }

  // Returns an invocable task that can be passed a cb and optionally a ctx
  // By defauult cb is a noop and ctx is a object literal
  return (cb = _.noop, ctx = {}) => {
    let taskState = lifecycleMethods.create(props)
    let operation

    if (lifecycleMethods.resumable(taskState)) {
      afterOperation(taskState)
    } else {
      taskState = lifecycleMethods.resolveInput(taskState)
      operation = lifecycleMethods.createOperation(taskState, operationCreator)

      lifecycleMethods.settleOperation(operation, (err, results) => {
        if (err) return afterOperation(taskState, err)

        return afterOperation(taskState, null, results)
      })
    }

    function afterOperation (taskState, err, results) {
      taskState = lifecycleMethods.resolveOutput(taskState)
      taskState = lifecycleMethods.validateOutput(taskState)
      taskState = lifecycleMethods.addDAGNode(taskState)

      cb(null, taskState)
    }
  }
}

module.exports = task
