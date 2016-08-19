'use strict'

// Effects
// Sagas is a redux middleware - it will watch for dispatched actions (which
// have already affected the store) and let you 'catch' them with the 'take'
// effect - effects are functions that pause the generator, letting the
// middleware running it do work on the yielded effect.
// call(fn, ...args) -> run a function asynchdronously, get back its resolved value
// put(action) -> dispatch an action to the store
// take(pattern) -> wait for an action matching `pattern`
// select(selector) -> select a piece of the store
// fork(fn) -> runs another generator in the background - used for tasks that first take()
const { buffers, eventChannel, END } = require('redux-saga')
const { call, put, take, select, fork } = require('redux-saga/effects')
const { EventEmitter2 } = require('eventemitter2')
const Promise = require('bluebird')

const { tab } = require('../utils')

// Mostly promises
// Ran as async side effects through call() in the middle of start/end/fail sagas
const {
  checkResumable,
  resolveInput,
  createOperation,
  settleOperation,
  resolveOutput,
  validateOutput,
  postValidation
} = require('../lifecycle')

// Store selector:
// (uid) => (store) => store.tasks[uid]
// Pass into select() effect
const { selectTask } = require('./selectors.js')

// Actions and action types
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
  successValidatatingOutput,
  SUCCESS_VALIDATING_OUTPUT,
  appendToLog,
  APPEND_TO_LOG
} = require('../reducers/tasks.js')

const { addOutput } = require('../reducers/collection.js')

/**
 * The task lifecycle.
 *
 * A series of generator functions, each which
 * 1. Waits on an action. Usually something like SUCCESS_PROPHASE
 * 2. Dispatches the START_METAPHASE action.
 * 3. Runs the main side effect(s) that are associated with that lifecycle stage.
 * 4. Dispatches either a SUCCESS_METAPHASE or FAIL_METAPHASE or ERR_METAPHASE
 * 5. The last action is assumed to be take()'d by another saga: else, the lifecycle ends
 * 6. (There is a potential for infinite looping - or recursive pipelines...)
 *
 * These "lifecycle method" chunks can be further composed, forked, etc.
 *
 * This lifecycle can be described as the following state machine:
 *
 *                                     (if resumable) ↘︎
 * resolveInput -> (!resumable) ->      operationSaga -> resolveOutput -> validateOutput
 *           ↖︎ (if resumable and output cannot be resolved)    ⏎
 *
 * The setting of resumable = false lets us avoid an infinite loop when output
 * cannot be resolved.
 *
 * @receives a CREATE_TASK action
 * @forks resolveInputSaga
 *  @receives START_RESOLVE_INPUT
 *  @emits SUCCESS/FAIL/ERR_RESOLVE_INPUT
 * @forks operationSaga
 *  @receives SUCCESS_RESOLVE_INPUT
 *  @emits ...
 * @forks resolveOutputSaga
 * @forks validateOutputSage
 */
function* lifecycle (action) {
  // Take task uid,
  // operationCreator: the actual task we are running, called with computing arguments
  // taskResolve, taskReject: from the parent promise which dispatched the action
  const { uid, operationCreator, taskResolve, taskReject } = action
  const logEmitter = new EventEmitter2({ wildcard: true })
  const nestedLogger = (depth) => ({
    emit: (channel, content, tabLevel = 0) => logEmitter.emit(channel, content, depth + tabLevel),
    depth
  })

  // These will all launch now but each blocks until ready
  // Each one waits for the success action from the preceding one
  // This is the regular (resume off or first-run) lifecyle, in order
  yield fork(resolveInputSaga, uid)
  yield fork(operationSaga, uid, operationCreator)
  yield fork(resolveOutputSaga, uid, 'after')
  yield fork(validateOutputSaga, uid)
  yield fork(postValidationSaga, uid)

  // Watches a channel made of logEmitter events
  yield fork(loggerSaga, uid)

  logEmitter.emit('log', `Running task ${uid}`)

  // If resumable, jump into resolveOutputSaga first (and skip waiting for operation)
  if (yield call(checkResumable, yield select(selectTask(uid)))) {
    yield* resolveOutputSaga(uid, 'before')
  } else {
    yield put(startResolveInput(uid))
  }

  function logger ({ uid }) {
    return eventChannel((emitter) => {
      logEmitter.on('*', function(content, tabLevel = 0) {
        // console.log('logger event triggered: ', content)
        emitter({
          uid,
          channel: this.event,
          content: !tabLevel ? content : tab(tabLevel) + content
        })
      })

      logEmitter.on('CLOSE_LOG', () => emitter(END))

      // Subscriber must return an unsubscribe function
      return () => logEmitter.removeAllListeners()
    }, buffers.sliding(100))
  }

  function* loggerSaga (uid) {
    const logChan = yield call(logger, { uid })
    try {
      while (true) {
        const { uid, channel, content } = yield take(logChan)
        yield put(appendToLog({ uid, channel, content }))
        // const currentLog = yield(select(s => s.tasks[uid].log.currentLog))
        // WHERE LOGGING HAPPENS
        console.log(content)
      }
    } finally {
      console.log('logger terminated and unsubscribed')
    }
  }

  function* resolveInputSaga (uid) {
    yield take(START_RESOLVE_INPUT)
    try {
      const results = yield call(resolveInput, yield select(selectTask(uid)), nestedLogger(1))

      yield put(successResolveInput(uid, results ? results.resolvedInput : null))
    } catch (err) {
      console.log('error: ', err)
    }
  }

  function* operationSaga (uid, operationCreator) {
    yield take(SUCCESS_RESOLVE_INPUT)
    let ops
    try {
      yield put(startOperation(uid))
      ops = yield call(createOperation, yield select(selectTask(uid)), operationCreator, nestedLogger(1))
      const settled = yield call(settleOperation, ops.operation)
      // console.log('post settled: ', settled)
      yield put(successOperation(uid))
    } catch (err) {
      console.log('err2: ', err)
    }
  }

  function* resolveOutputSaga (uid, mode) {
    if (mode !== 'before') yield take(SUCCESS_OPERATION)

    yield put(startResolveOutput(uid))
    try {
      const resolvedOutput = yield call(resolveOutput, yield select(selectTask(uid)), nestedLogger(1))
      yield put(successResolveOutput(uid, resolvedOutput))
    } catch (err) {
      console.log('error3: ', err)

      if (mode === 'before') yield put(startResolveInput(uid))
    }
  }

  function* validateOutputSaga (uid) {
    yield take(SUCCESS_RESOLVE_OUTPUT)

    yield put(startValidatingOutput(uid))
    try {
      const validations = yield call(validateOutput, yield select(selectTask(uid)), nestedLogger(1))
      yield put(successValidatatingOutput(uid))
    } catch (err) {
      console.error('err4: ', err)
    }

    yield call(() => logEmitter.emit('log', `Finished task ${uid}`))
    // logEmitter.emit('CLOSE_LOG')
    yield put(successResolveOutput(uid))
  }

  function* postValidationSaga (uid) {
    yield take(SUCCESS_VALIDATING_OUTPUT)
    logEmitter.emit('log', tab(1) + 'Adding task output to DAG')
    yield put(addOutput(uid, yield select(selectTask(uid))))

    yield* finish(uid)
  }

  function* finish (uid) {
    // Resolve the promise from task.js
    taskResolve(yield select(selectTask(uid)))
  }
}


module.exports = lifecycle
