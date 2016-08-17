'use strict'

const { call, put, take, select, fork } = require('redux-saga/effects')

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

const {
  startResolveInput,
  START_RESOLVE_INPUT,
  successResolveInput,
  SUCCESS_RESOLVE_INPUT,
  startOperation,
  successOperation,
  SUCCESS_OPERATION,
  startResolveOutput,
  successResolveOutput,
  SUCCESS_RESOLVE_OUTPUT,
  startValidatingOutput,
  successValidatatingOutput
} = require('../reducers/tasks.js')

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
  const { uid, operationCreator, taskResolve, taskReject } = action

  // These will all launch now but each blocks until ready
  // Each one waits for the success action from the preceding one
  // This is the regular (resume off or first-run) lifecyle, in order
  yield fork(resolveInputSaga, uid)
  yield fork(operationSaga, uid, operationCreator)
  yield fork(resolveOutputSaga, uid, 'after')
  yield fork(validateOutputSaga, uid)

  // If resumable, jump into resolveOutputSaga first (and skip waiting for operation)
  if (yield call(checkResumable, yield select(selectTask(uid)))) {
    yield* resolveOutputSaga(uid, 'before')
  } else {
    yield put(startResolveInput(uid))
  }

  function* resolveInputSaga (uid) {
    yield take(START_RESOLVE_INPUT)
    console.log('resolve input triggered')
    try {
      const resolvedInput = yield call(resolveInput, yield select(selectTask(uid)))
      console.log('_saga: resolvedInput: ', resolvedInput)
      yield put(successResolveInput(uid, resolvedInput))
      // console.log(yield select(selectTask(uid)))
    } catch (err) {
      console.log('error: ', err)
    }
  }

  function* operationSaga (uid, operationCreator) {
    yield take(SUCCESS_RESOLVE_INPUT)
    console.log('operation saga triggered')
    let ops
    try {
      yield put(startOperation(uid))
      ops = yield call(createOperation, yield select(selectTask(uid)), operationCreator)
      const settled = yield call(settleOperation, ops.operation)
      console.log('post settled: ', settled)
      yield put(successOperation(uid))
    } catch (err) {
      console.log('err2: ', err)
    }
  }

  function* resolveOutputSaga (uid, mode) {
    if (mode !== 'before') yield take(SUCCESS_OPERATION)

    console.log('resolve output triggered')

    yield put(startResolveOutput(uid))
    try {
      const resolvedOutput = yield call(resolveOutput, yield select(selectTask(uid)))
      yield put(successResolveOutput(uid, resolvedOutput))
    } catch (err) {
      console.log('error3: ', err)

      if (mode === 'before') yield put(startResolveInput(uid))
    }
  }

  function* validateOutputSaga (uid) {
    yield take(SUCCESS_RESOLVE_OUTPUT)
    console.log('validate triggered')

    yield put(startValidatingOutput(uid))
    try {
      const validations = yield call(validateOutput, yield select(selectTask(uid)))
      console.log('validations: ', validations)
      yield put(successValidatatingOutput(uid))
    } catch (err) {
      console.error('err4: ', err)
    }

    taskResolve(yield select(selectTask(uid)))
  }
}


module.exports = lifecycle
