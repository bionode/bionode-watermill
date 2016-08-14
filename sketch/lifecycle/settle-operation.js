'use strict'

const _ = require('lodash')
const asyncDone = require('async-done')

const settleOperation = (operation, cb = _.noop) => {
  asyncDone(
    () => operation,
    cb
  )
}


module.exports = settleOperation
