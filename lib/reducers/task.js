'use strict'

const stringify = require('json-stable-stringify')
const defaultConfig = require('./config.js').defaultState
const { hash } = require('../utils.js')

// Actions
const CREATE_TASK = 'task/create'

// Statuses
const STATUS_CREATED = 'STATUS_CREATED'

// null values NEED to be set
// except for container, resolvedInput, resolvedOuput
// TODO a "deepObjectAssign", so that properties like params can be extended
const defaultState = {
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

const reducer = (state = defaultState, action) => {
  switch(action.type) {
    case CREATE_TASK:
      const task = Object.assign({}, state, action.task, { created: Date.now() })
      const hashes = {
        input: hash(stringify(task.input)),
        output: hash(stringify(task.output)),
        params: hash(stringify(task.params))
      }
      const uid = hash(hashes.input + hashes.output + hashes.params)
      Object.assign(task, { hashes, uid })
      return task
      break
    default:
      return state
  }
}

reducer._createTask = (task) => reducer(undefined, {
  type: CREATE_TASK,
  task
})

module.exports = reducer
