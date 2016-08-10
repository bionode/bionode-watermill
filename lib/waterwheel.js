'use strict'

const { createStore, applyMiddleware } = require('redux')
const thunk = require('redux-thunk').default
const logger = require('./utils/my-redux-logger.js')
const rootReducer = require('./reducers')

const { shell } = require('./utils/shell.js')
const parallel = require('./orchestration/parallel.js')
const join = require('./orchestration/join.js')

// Make a waterwheel instance
const waterwheel = () => {
  const store = createStore(rootReducer, applyMiddleware(thunk/*, logger*/))

  const actions = require('./utils/get-actions.js')(store.dispatch)

  const task = require('./orchestration/task-runner.js')(store, actions)

  return {
    task,
    join,
    parallel,
    shell,
    store
  }
}

module.exports = waterwheel

