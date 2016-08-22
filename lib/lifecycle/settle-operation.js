'use strict'

const _ = require('lodash')
const asyncDone = require('async-done')

const settleOperation = (operation) => new Promise((resolve, reject) => {
  asyncDone(
    () => operation,
    (err, data) => err ? reject(err) : resolve(data)
  )
})


module.exports = settleOperation
