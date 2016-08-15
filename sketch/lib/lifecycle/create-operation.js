'use strict'

const _ = require('lodash')
const path = require('path')
const Promise = require('bluebird')
const mkdirp = Promise.promisify(require('mkdirp'))
const fs = Promise.promisifyAll(require('fs'))

const { shell } = require('../utils/shell.js')

/**
 * Produce operation parameters give the current task state.
 *
 * @param props {Object}
 */
const createOperationProps = (taskState) => new Promise((resolve, reject) => {
  const dir = taskState.dir

  if (taskState.resolvedInput) {
    const symlinkTarget = taskState.resolvedInput
    const symlinkPath = path.resolve(dir, path.basename(symlinkTarget))

    fs.symlinkAsync(symlinkTarget, symlinkPath).then(() => {
      const operationProps = {
        input: symlinkPath,
        output: taskState.resolvedOutput
      }
      resolve(Object.assign({}, taskState, operationProps))
    }).catch(err => {
      console.log('got symlink err: ', err)
    })
  } else {
    const operationProps = {
      input: taskState.resolvedInput,
      output: taskState.resolvedOutput
    }
    resolve(Object.assign({}, taskState, operationProps))
  }
})

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
    createOperationProps(taskState).then((operationProps) => {
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
})

module.exports = createOperation
