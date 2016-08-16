'use strict'

const { takeEvery } = require('redux-saga')
const { call, put, take, select } = require('redux-saga/effects')

/**
 * The task lifecycle
 *
 * `call`s all side effects and `put`s some actions.
 * This lets us seperate the async side effects from their results and
 * application of those results to the store, while also providing an iterable
 * that is easy to test.
 */
function* lifecycle (action) {
  const { operationCreator, taskResolve, taskReject } = action
  console.log('state from sagas: ', JSON.stringify(yield select(), null, 2))

  console.log('operationCreator: ', action.operationCreator)

  console.log('gonna resolve you son')
  taskResolve('foo')
}

function* rootSaga() {
  yield* takeEvery('tasks/create', lifecycle)
}

// TODO differnet modules
rootSaga.lifecycle = lifecycle

module.exports = rootSaga

