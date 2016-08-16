'use strict'

const { call, put, take, select } = require('redux-saga/effects')

const {
  checkResumable
} = require('../lifecycle')
const { selectTask } = require('./selectors.js')

/**
 * The task lifecycle
 *
 * `call`s all side effects and `put`s some actions.
 * This lets us seperate the async side effects from their results and
 * application of those results to the store, while also providing an iterable
 * that is easy to test.
 *
 * Receives a CREATE_TASK
 */
function* lifecycle (action) {
  const { operationCreator, taskResolve, taskReject } = action

  const taskState = yield select(selectTask(action.uid))

  const resumable = yield call(checkResumable, taskState)

  console.log('operationCreator: ', action.operationCreator)

  console.log('gonna resolve you son')
  taskResolve('foo')
}

module.exports = lifecycle
