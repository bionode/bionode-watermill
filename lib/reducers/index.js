'use strict'

const { combineReducers } = require('redux')

const config = require('./config.js')
const tasks = require('./tasks.js')

module.exports = combineReducers({
  config,
  tasks
})
