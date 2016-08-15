'use strict'

const _ = require('lodash')
const path = require('path')
const Promise = require('bluebird')
const mkdirp = Promise.promisify(require('mkdirp'))

const { shell } = require('../utils/shell.js')

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
  return Object.assign({}, taskState, operationProps)
}

/**
 * Create operation. Takes taskState.operationProps
 */
const createOperation = (taskState, operationCreator) => new Promise((resolve, reject) => {
  console.log('resolvedInput before making: ', taskState.resolvedInput)
  const dir = taskState.dir
  const madeDir = mkdirp(dir)
    .then(() => console.log('Created ' + dir))
    .catch(err => reject(err))

  madeDir.then(() => {
    const operationProps = createOperationProps(taskState)
    const operation = operationCreator(operationProps)

    // Convert string or array of strings into shell
    if (_.isString(operation)) {
      resolve({
        operation: shell(operation, { cwd: dir }),
        dir,
        operationString: operation
      })
    } else {
      resolve({
        operation,
        dir,
        operationString: JSON.stringify(operation)
      })
    }
  })
})

module.exports = createOperation
