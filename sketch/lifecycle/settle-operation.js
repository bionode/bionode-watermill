'use strict'

const _ = require('lodash')
const asyncDone = require('async-done')

const settleOperation = (operation, cb) => {
  if (!_.isFunction(cb)) throw new Error('No callback provided to settleOperation(operation, cb)')

  asyncDone(
    operation,
    cb
  )
}


module.exports = settleOperation
