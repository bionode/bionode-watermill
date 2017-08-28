'use strict'

const path = require('path')

// Reducers - this implies there is reducer proxying in this reducer
// const taskReducer = require('./task.js')
const { defaultTask } = require('../constants/default-task-state.js')
const config = require('../constants/default-config-state.js')

const { tab } = require('../utils')

// Actions
const CREATE_TASK = 'tasks/create'
const START_RESOLVE_INPUT = 'task/start_resolve-input'
const SUCCESS_RESOLVE_INPUT = 'task/success_resolve-input'

const START_OPERATION = 'task/start_operation'
const SUCCESS_OPERATION = 'task/success_operation'

const START_RESOLVE_OUTPUT = 'task/start_resolve-output'
const SUCCESS_RESOLVE_OUTPUT = 'task/success_resolve-output'

const START_VALIDATING_OUTPUT = 'task/start_validating-output'
const SUCCESS_VALIDATING_OUTPUT = 'task/success_validating-output'

const APPEND_TO_LOG = 'task/append-to-log'

const SET_OPERATION_STRING = 'task/set-operation-string'

// Internal Methods

/**
 * updateTask
 * Return a new state with a new task. Update a key in a dictionary (object)
 * with a new value, returning a new, referentially dis-equal object.
 * @param state {Object} reference to current state
 * @param task {Object} the task to add/update into state. Must have `task.uid`.
 */
const updateTask = (state, task) => {
  const { uid } = task

  if (!uid) {
    throw new Error('Cannot updateTask without a uid')
    return state
  }

  return Object.assign({}, state, { [uid]:task })
}

/**
 * @param status {string}
 * @returns {Number}
 */
const statusToTabLevel = (status) => {
  switch (status) {
    case 'RESOLVING_INPUT':
      return 1
      break
    default:
      return 0
  }
}

// Is never triggered with an undefined action:
// Called as a nested reducer from the tasks dictionary
const taskReducer = (state = {}, action) => {
  switch (action.type) {
    case START_RESOLVE_INPUT:
      return Object.assign({}, state, { status: 'RESOLVING_INPUT' })
      break
    case SUCCESS_RESOLVE_INPUT:
      return Object.assign({}, state, { resolvedInput: action.resolvedInput })
    case START_OPERATION:
      return Object.assign({}, state, { status: 'RUNNING_OPERATION' })
      break
    case START_RESOLVE_OUTPUT:
      return Object.assign({}, state, { status: 'RESOLVING_OUTPUT' })
      break
    case SUCCESS_RESOLVE_OUTPUT:
      return Object.assign({}, state, { resolvedOutput: action.resolvedOutput })
    case START_VALIDATING_OUTPUT:
      return Object.assign({}, state, { status: 'VALIDATING_OUTPUT' })
      break
    case SUCCESS_VALIDATING_OUTPUT:
      return Object.assign({}, state, { status: 'POST_VALIDATION', validated: true })
      break
    case SET_OPERATION_STRING:
      return Object.assign({}, state, { operationString: action.operationString })
      break
    case APPEND_TO_LOG:
      const { channel, content } = action
      const currentLog = Buffer.from(content)
      // console.log('currentLog:' + currentLog)
      // console.log(currentLog.toString())
      const newBuf = Buffer.concat([
        state.log[channel],
        Buffer.from('\n'),
        currentLog
      ], state.log[channel].length + 1 + currentLog.length)
      const newLog = Object.assign({}, state.log, { [channel]:newBuf, currentLog })
      return Object.assign({}, state, { log: newLog })
      break
    default:
      return state
  }
}

// The `tasks` store item is an object keyed by `task` uids. Each uid points to
// a `task` instance, which also contains `task.uid` for  redundancy. Thus, the
// default state is no tasks.
const defaultState = {}

const reducer = (state = defaultState, action) => {
  const { uid, type } = action

  if (type.startsWith('task/')) {
    // An action has been dispatched for an *assumed existing* task.
    // Nest the task reducer behind this one.

    // Get the action uid so we can retrieve the current state of that task.
    // If it is not provided, throw and return current state.
    const { uid } = action
    if (!uid) {
      throw new Error('task action not provided with uid')
      return state
    }

    // Use uid to retrieve current task. Throw and return if it is not found.
    const currentTask = state[uid]
    if (!currentTask) {
      throw new Error('task uid does not match an existing task')
      return state
    }

    // Call taskReducer with current state of corresponding task.
    // This is the point of proxying/nesting.
    const newTask = taskReducer(currentTask, action)

    return updateTask(state, newTask)
  }

  switch(type) {
    case CREATE_TASK:
      // Double check values being passed from action creator
      const { hashes, props, context } = action
      const newTask = Object.assign(
        {},
        defaultTask,
        { dir: path.resolve(defaultTask.dir, uid.substring(0, 7)) },
        // User can overwrite dir
        props,
        { uid, hashes },
        { created: Date.now() }
      )
      newTask.context = context
      const newState = updateTask(state, newTask)

      return newState
    default:
      return state
  }
}

reducer.CREATE_TASK = CREATE_TASK // exported constants along functions is nicer with JSM..
reducer.createTask = ({
  uid, task, props, hashes, context, operationCreator, taskResolve, taskReject
}) => ({
  type: CREATE_TASK,
  uid, task, props, hashes, context, operationCreator, taskResolve, taskReject
})

reducer.startResolveInput = (uid) => ({ type: START_RESOLVE_INPUT, uid })
reducer.successResolveInput = (uid, resolvedInput) => ({ type: SUCCESS_RESOLVE_INPUT, uid, resolvedInput })
reducer.SUCCESS_RESOLVE_INPUT = SUCCESS_RESOLVE_INPUT

reducer.startOperation = (uid) => ({ type: START_OPERATION, uid })
reducer.successOperation = (uid) => ({ type: SUCCESS_OPERATION, uid })
reducer.SUCCESS_OPERATION = SUCCESS_OPERATION

// passes operationString to action and next taskState
reducer.setOperationString = (uid, operationString) => ({ type:SET_OPERATION_STRING, uid, operationString })
reducer.SET_OPERATION_STRING = SET_OPERATION_STRING

reducer.startResolveOutput = (uid) => ({ type: START_RESOLVE_OUTPUT, uid })
reducer.START_RESOLVE_INPUT = START_RESOLVE_INPUT
reducer.successResolveOutput = (uid, resolvedOutput) => ({ type: SUCCESS_RESOLVE_OUTPUT, uid, resolvedOutput })
reducer.SUCCESS_RESOLVE_OUTPUT = SUCCESS_RESOLVE_OUTPUT

reducer.startValidatingOutput = (uid) => ({ type: START_VALIDATING_OUTPUT, uid })
reducer.successValidatatingOutput = (uid) => ({ type: SUCCESS_VALIDATING_OUTPUT, uid })
reducer.SUCCESS_VALIDATING_OUTPUT = SUCCESS_VALIDATING_OUTPUT


// reducer.appendToLog = ({ uid, channel = 'log', content }) => {
//   const tabLevel =
// }

reducer.appendToLog = ({ uid, channel = 'log', content }) => ({
  type: APPEND_TO_LOG,
  uid, channel, content
})
reducer.APPEND_TO_LOG = APPEND_TO_LOG

module.exports = reducer
