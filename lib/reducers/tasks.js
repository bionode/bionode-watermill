'use strict'

const stringify = require('json-stable-stringify')
const defaultConfig = require('./config.js').defaultState
const { hash } = require('../utils.js')

// Actions
const CREATE_TASK = 'tasks/create'
const CHECK_RESUMABLE = 'tasks/check-resumable'
const RESOLVE_OUTPUT = 'tasks/resolve-output'

// Statuses
const {
  STATUS_CREATED,
  CHECKING_RESUMABLE,
  RESOLVING_INPUT
} = require('../constants/task-status-types.js')

const initTask = (task) => {
  // null values NEED to be set
  // except for container, resolvedInput, resolvedOuput
  // TODO a "deepObjectAssign", so that properties like params can be extended
  const defaultTask = {
    threads: defaultConfig.threads,
    container: defaultConfig.container,
    resume: defaultConfig.resume,
    uid: null,
    hashes: {
      input: null,
      output: null,
      params: null
    },
    name: 'Unnamed Task',
    dir: null,
    input: null,
    output: null,
    status: STATUS_CREATED,
    created: null,
    params: {}
  }

  task = Object.assign({}, defaultTask, task, { created: Date.now() })
  const hashes = {
    input: hash(stringify(task.input)),
    output: hash(stringify(task.output)),
    params: hash(stringify(task.params))
  }
  const uid = hash(hashes.input + hashes.output + hashes.params)
  Object.assign(task, { hashes, uid })
  return task
}

const defaultState = {}
const reducer = (state = defaultState, action) => {
  switch (action.type) {
    case CREATE_TASK:
      return createTaskHandler(state, action)
      break
    case CHECK_RESUMABLE:
      return checkResumableHandler(state, action)
      break
    case RESOLVE_OUTPUT:
      return resolveOutputHandler(state, action)
      break
    default:
      return state
  }
}

function createTaskHandler(state, { task }) {
  return Object.assign({}, state, { [task.uid]:task })
}

reducer.createTask = (task) => (dispatch) => new Promise((resolve, reject) => {
  const newTask = initTask(task)

  dispatch({
    type: CREATE_TASK,
    task: newTask
  })

  resolve(newTask.uid)
})

function checkResumableHandler(state, { uid }) {
  const task = state[uid]
  const { resume } = task

  if (resume === 'off') {
    // No resuming
    task.status = RESOLVING_INPUT
  } else if (resume === 'on') {
    // Going to initiate resume attempt
    task.status = CHECKING_RESUMABLE
  }

  return Object.assign({}, state, { [uid]:task })
}

reducer.checkResumable = (uid) => (dispatch) => new Promise((resolve, reject) => {
  dispatch({
    type: CHECK_RESUMABLE,
    uid
  })

  resolve(uid)
})

function resolveOutputHandler(state, { uid, resolvedOutput }) {
  const task = state[uid]

  task.resolvedOutput = resolvedOutput

  return Object.assign({}, state, { [uid]:task })
}

reducer.resolveOutput = (task) => (dispatch) => new Promise((resolve, reject) => {
  // TODO more actions
  const outerResolver = require('../resolvers.js')
  const validators = require('../validators.js')

  const { output } = task

  let alreadyOutput = outerResolver(output, task, 'checking', {
    fileResolvers: [validators.matchHash, validators.matchTime]
  })

  // Check that everything resolved
  let safe = true
  if (!(alreadyOutput instanceof Array)) {
    alreadyOutput = [alreadyOutput]
  }

  for (let item of alreadyOutput) {
    if (!item.resolved) {
      safe = false
    }
  }

  if (safe) {
    // TODO handle logging better
    console.log('  ' + 'Skipping this task')
    console.log('todo: kill stream after setting its ouput')

    dispatch({
      type: RESOLVE_OUTPUT,
      uid: task.uid,
      resolvedOutput: alreadyOutput,
    })

    // task.stream.output = () => {
    //   return alreadyOutput
    // }
    //
    // task.stream.destroy()
  } else {
    console.log('TODO not safe')
  }

  resolve(task.uid)
})

module.exports = reducer
