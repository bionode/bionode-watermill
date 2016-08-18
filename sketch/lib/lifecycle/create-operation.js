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
const createOperation = (taskState, operationCreator, logger) => new Promise((resolve, reject) => {
  logger.emit('log', 'Creating operation')
  // console.log('resolvedInput before making: ', taskState.resolvedInput)
  const dir = taskState.dir
  const madeDir = mkdirp(dir)
    .then(() => logger.emit('log', 'Created ' + dir))
    .catch(err => reject(err))

  madeDir.then(() => {
    createOperationProps(taskState).then((operationProps) => {
      const operation = operationCreator(operationProps)

      logger.emit('log', 'operation as string: ' + operationCreator.toString())

      // Convert string or array of strings into shell
      if (_.isString(operation)) {
        resolve({
          operation: shell(operation, { cwd: dir }, {
            // Wrap logger to introduce default tabLevel
            emit: (channel, content, tabLevel = 0) => logger.emit(channel, content, logger.depth + tabLevel)
          }),
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
