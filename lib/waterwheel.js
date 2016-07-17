'use strict'

const { createStore, applyMiddleware } = require('redux')
const thunk = require('redux-thunk').default
const logger = require('./my-redux-logger.js')
const rootReducer = require('./reducers')

const store = createStore(rootReducer, applyMiddleware(thunk, logger))

const actions = require('./get-actions.js')(store.dispatch)

const task = require('./Task.js')(actions)
const join = require('./Join.js')
const parallel = require('./parallel.js')
const { shell } = require('./wrappers.js')

module.exports = {
  task,
  join,
  parallel,
  shell,
  store
}
