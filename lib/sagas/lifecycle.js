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
  setupTaskdir,
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
  APPEND_TO_LOG,
  setOperationString,
  SET_OPERATION_STRING
} = require('../reducers/tasks.js')

const { addOutput } = require('../reducers/collection.js')

/**
 * The task lifecycle.
 *
 * TODO rewrite this as per now using yield*
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
 */
function* lifecycle (action) {
  // Take task uid,
  // operationCreator: the actual task we are running, called with computing arguments
  // taskResolve, taskReject: from the parent promise which dispatched the action
  // taskResolve will trigger the myTask().then, similary taskReject -> catch
  const { uid, operationCreator, taskResolve, taskReject } = action
  const miniUid = uid.substring(0, 7)

  // Logging - pass an emitter into functions that can log for the task
  // This is more verbose than console.log - but lets us manage interpolated
  // logging for different tasks running at once. Can either print logs
  // interpolated, prefixed with task uid, or after completion of each task
  // Alternatively, could pass store.dispatch into each function, but using a
  // saga lets us manage all logs in one place.
  const logEmitter = new EventEmitter2({ wildcard: true })
  const nestedLogger = (depth) => ({
    emit: (channel, content, tabLevel = 0) => logEmitter.emit(channel, content, depth + tabLevel),
    depth
  })
  // Watches a channel made of logEmitter events
  yield fork(loggerSaga, uid)

  // This task will not have resolveInput, resolveOutput (and never will)
  const originalTask = yield select(selectTask(uid))
  logEmitter.emit('log', `=== ${originalTask.name} ===`)
  logEmitter.emit('log', `Running task ${miniUid}`)

  function* orderedLifecycle () {
    let next
    next = yield* resolveInputSaga()
    if (next) next = yield* operationSaga(operationCreator)
    if (next) next = yield* afterResolveOutput({})
  }

  function* afterResolveOutput ({ skipResolve = false }) {
    let next = true
    // Skip resolve for resumables that already resolved output
    if (!skipResolve) next = yield* resolveOutputSaga()
    if (next) next = yield* validateOutputSaga()
    if (next) next = yield* postValidationSaga()
    if (next) {
      const { delay } = require('redux-saga')
      // Wait a bit for logs to catch up
      // TODO fix
      yield delay(500)
      yield* finish(action.uid)
    }
  }

  if (originalTask.resume === 'off') {
    yield* orderedLifecycle()
  } else if (originalTask.resume === 'on') {
    // Will trigger orderedLifecycle if mode === 'before' on failure to resolveOutput
    yield* resolveOutputSaga('before')
  } else {
    logEmitter.log('log', 'WARN: resume was not "on" or "off"')
    yield* orderedLifecycle()
  }

  function* resolveInputSaga () {
    try {
      const results = yield call(
        resolveInput,
        yield select(selectTask(uid)),
        yield select(s => s.collection),
        nestedLogger(1)
      )
      yield put(successResolveInput(uid, results ? results.resolvedInput : null))
      return true
    } catch (err) {
      console.log('error!: ', err.toString())

      yield* fail(uid)
      return false
    }
  }

  function* operationSaga (operationCreator) {
    let ops
    try {
      yield put(startOperation(uid))
      const linkedResolvedInput = yield call(setupTaskdir, yield(select(selectTask(uid))), nestedLogger(1))
      ops = yield call(
        createOperation,
        yield select(selectTask(uid)), linkedResolvedInput, operationCreator, nestedLogger(1)
      )
      // console.log('operation string: ', ops.operationString)
      // this pushes operationString to taskState
      yield put(setOperationString(uid, ops.operationString))
      const settled = yield call(settleOperation, ops.operation)
      // console.log('post settled: ', settled)
      yield put(successOperation(uid))
      return true
    } catch (err) {
      console.log('err2: ', err.toString())

      yield* fail(uid)
      return false
    }
  }

  function* resolveOutputSaga (mode) {
    yield put(startResolveOutput(uid))

    try {
      const resolvedOutput = yield call(resolveOutput, yield select(selectTask(uid)), nestedLogger(1))
      logEmitter.emit('log', 'resolvedOutput: ' + resolvedOutput)
      yield put(successResolveOutput(uid, resolvedOutput))

      if (mode === 'before') yield* afterResolveOutput({ skipResolve: true })

      return true
    } catch (err) {
      logEmitter.emit('log', err.toString())

      if (mode === 'before') yield* orderedLifecycle()

      yield* fail(uid)
      return false
    }
  }

  function* validateOutputSaga () {
    yield put(startValidatingOutput(uid))

    try {
      const validations = yield call(validateOutput, yield select(selectTask(uid)), nestedLogger(1))
      yield put(successValidatatingOutput(uid))

      return true
    } catch (err) {
      console.error('err4: ', err.toString())

      yield* fail(uid)
      return false
    }
  }

  function* postValidationSaga () {
    logEmitter.emit('log', tab(1) + `Adding task output to DAG for ${uid}`)
    try {
      yield put(addOutput(uid, yield select(selectTask(uid))))
      const DAG = yield(select(s => s.collection))
      const { jsonifyGraph } = require('../reducers/collection.js')
      logEmitter.emit('log', 'current DAG: ' + JSON.stringify(jsonifyGraph(DAG), null, 2))

      return true
    } catch (err) {
      logEmitter.emit('log', err.toString())

      yield* fail(uid)
      return false
    }
  }

  function* finish (uid) {
    // Resolve the promise from task.js
    const finalTask = yield select(selectTask(uid))
    // TODO if was ran alongside other tasks
    // console.log(task.log.log.toString())
    logEmitter.emit('log', `Finished task ${miniUid}`)
    taskResolve(finalTask)
  }

  function* fail (uid) {
    const finalTask = yield select(selectTask(uid))
    taskReject(finalTask)
  }

  function logger () {
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
    const logChan = yield call(logger)
    try {
      while (true) {
        const { uid, channel, content } = yield take(logChan)
        yield put(appendToLog({ uid, channel, content }))
        // const currentLog = yield(select(s => s.tasks[uid].log.currentLog))
        // WHERE LOGGING HAPPENS
        // TODO trigger based on join, parallel, fork
        console.log(miniUid, ':', content)
      }
    } finally {
      console.log('Logger terminated and unsubscribed for task ' + miniUid)
    }
  }
}


module.exports = lifecycle
