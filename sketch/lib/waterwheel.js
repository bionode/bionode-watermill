'use strict'

const task = require('./task.js')
const join = require('./orchestrators/join.js')
const parallel = require('./orchestrators/parallel.js')

module.exports = {
  task,
  join,
  parallel
}
