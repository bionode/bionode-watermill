'use strict'

const { combineReducers } = require('redux')

const config = require('./config.js')
const tasks = require('./tasks.js')
const task = require('./task.js')

module.exports = combineReducers({
  config,
  tasks,
  task
})
