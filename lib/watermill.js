'use strict'

const { createStore, applyMiddleware } = require('redux')
const createLogger = require('redux-logger')
const thunk = require('redux-thunk').default
const rootReducer = require('./reducers')

const createSagaMW = require('redux-saga').default
const sagaMW = createSagaMW()
const logger = createLogger()

// default middlewares
let middlewares = [thunk, sagaMW]

// if Environmental variable REDUX_LOGGER is 1 then...
const loggerEnv = process.env.REDUX_LOGGER
// for now call REDUX_LOGGER=1 before executing a pipeline
// REDUX_LOGGER=1 node pipeline.js
if (loggerEnv === '1') {
  middlewares.push(logger)
}

const store = createStore(rootReducer, applyMiddleware(...middlewares))

const rootSaga = require('./sagas')
sagaMW.run(rootSaga)

const task = require('./task.js')(store.dispatch)
const join = require('./orchestrators/join.js')(store.dispatch)
const junction = require('./orchestrators/junction.js')(store.dispatch)
const fork = require('./orchestrators/fork.js')

module.exports = {
  task,
  join,
  junction,
  fork,
  store
}
