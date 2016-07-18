'use strict'

const stringify = require('json-stable-stringify')
const defaultConfig = require('./config.js').defaultState
const { hash } = require('../utils.js')

// Actions
const CREATE_TASK = 'tasks/create'

// Statuses
const STATUS_CREATED = 'STATUS_CREATED'

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
      const { task } = action
      console.log(task.uid)
      return Object.assign({}, state, { [task.uid]:task })
      break
    default:
      return state
  }
}

reducer.createTask = (task) => (dispatch) => new Promise((resolve, reject) => {
  const newTask = initTask(task)

  dispatch({
    type: CREATE_TASK,
    task: newTask
  })

  resolve(newTask.uid)
})

module.exports = reducer
