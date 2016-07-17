'use strict'

const task = require('./Task.js')
const join = require('./Join.js')
const parallel = require('./parallel.js')
const { shell } = require('./wrappers.js')

const { createStore, applyMiddleware } = require('redux')
const thunk = require('redux-thunk').default
const logger = require('./my-redux-logger.js')

const rootReducer = require('./reducers')

const store = createStore(rootReducer, applyMiddleware(thunk, logger))

module.exports = {
  task,
  join,
  parallel,
  shell,
  store
}
