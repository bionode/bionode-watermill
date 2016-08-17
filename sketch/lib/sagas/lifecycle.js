'use strict'

const { call, put, take, select } = require('redux-saga/effects')

const {
  checkResumable,
  resolveInput,
  createOperation,
  settleOperation,
  resolveOutput,
  validateOutput,
  postValidation
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
  if (resumable) {
    yield* afterOperation(action.uid)
  }

  // console.log('operationCreator: ', action.operationCreator)
  // console.log('gonna resolve you son')
  // taskResolve('foo')

  function* beforeOperation (uid) {
    const resolvedInput = yield call(resolveInput, yield select(selectTask(uid)))
    const { operation } = yield call(
      createOperation,
      Object.assign(yield select(selectTask(uid))),
      operationCreator
    )
    yield call(settleOperation, operation)

    yield* afterOperation(uid)
  }

  function* afterOperation (uid) {
    let resolvedOutput
    try {
      resolvedOutput = yield call(resolveOutput, yield select(selectTask(uid)))
      console.log('resolvedOutput: ', resolvedOutput)

      yield call(
        validateOutput,
        Object.assign(yield select(selectTask(uid)), { resolvedOutput: resolvedOutput.resolvedOutput })
      )
      taskResolve('resolve from saga')
    } catch(err) {
      console.log('err: ', err)

      // TODO fix this inf. loop that can happen here
      yield* beforeOperation(uid)
    }
  }
}


module.exports = lifecycle
