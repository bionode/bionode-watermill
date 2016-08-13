'use strict'

const asyncDone = require('async-done')

/**
 * Create first taskState given props
 */
const create = (props) => {
  return Object.assign({}, props, { created: Date.now() })
}

const resolveInput = (taskState) => Object.assign({}, taskState)

const resumable = (taskState) => taskState.resumbable === 'on' ? true : false

/**
 * Produce operation parameters give the current task state.
 *
 * @param props {Object}
 */
const createOperationProps = (taskState) => {
  const operationProps = {
    input: taskState.resolvedInput,
    output: taskState.resolvedOutput
  }
  return Object.assign({}, taskState, { operationProps })
}

/**
 * Create operation. Takes taskState.operationProps
 */
const createOperation = (taskState, operationCreator) => {
  const operationProps = createOperationProps(taskState)

  return operationCreator(operationProps)
}

const settleOperation = (operation, cb) => asyncDone(() => operation, cb)

const resolveOutput = (taskState) => Object.assign({}, taskState)

const validateOutput = (taskState) => Object.assign({}, taskState)

const addDAGNode = (taskState) => Object.assign({}, taskState)

const lifecycleMethods = {
  create,
  resumable,
  resolveInput,
  createOperation,
  settleOperation,
  resolveOutput,
  validateOutput,
  addDAGNode
}

module.exports = lifecycleMethods
