'use strict'

const { createStore, applyMiddleware } = require('redux')
const thunk = require('redux-thunk').default
const rootReducer = require('./reducers')

const store = createStore(rootReducer, applyMiddleware(thunk))

const task = require('./task.js')(store)
const join = require('./orchestrators/join.js')
const parallel = require('./orchestrators/parallel.js')

module.exports = {
  task,
  join,
  parallel,
  store
}
