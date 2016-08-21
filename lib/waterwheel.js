'use strict'

const { createStore, applyMiddleware } = require('redux')
const thunk = require('redux-thunk').default
const rootReducer = require('./reducers')

const createSagaMW = require('redux-saga').default
const sagaMW = createSagaMW()

const store = createStore(rootReducer, applyMiddleware(thunk, sagaMW))

const rootSaga = require('./sagas')
sagaMW.run(rootSaga)

const task = require('./task.js')(store.dispatch)
const join = require('./orchestrators/join.js')
const junction = require('./orchestrators/junction.js')(store.dispatch)

module.exports = {
  task,
  join,
  junction,
  store
}
