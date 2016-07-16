'use strict'

const task = require('./Task.js')
const join = require('./Join.js')
const parallel = require('./parallel.js')
const { shell } = require('./wrappers.js')

module.exports = {
  task,
  join,
  parallel,
  shell
}
