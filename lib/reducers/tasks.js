'use strict'

const stringify = require('json-stable-stringify')
const defaultConfig = require('./config.js').defaultState
const { hash } = require('../utils.js')

// Actions
const CREATE_TASK = 'tasks/create'
const CHECK_RESUMABLE = 'tasks/check-resumable'
const RESOLVE_OUTPUT = 'tasks/resolve-output'
const RESOLVE_INPUT = 'tasks/resolve-input'
const CREATE_ACTION = 'tasks/create-action'

// Statuses
const {
  STATUS_CREATED,
  CHECKING_RESUMABLE,
  RESOLVING_INPUT
} = require('../constants/task-status-types.js')

const initTask = (task) => {
  // null values NEED to be set
  // except for container, resolvedInput, resolvedOuput, action
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
    action: null,
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
    case RESOLVE_INPUT:
      return resolveInputHandler(state, action)
    case CREATE_ACTION:
      return createActionHandler(state, action)
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
  } else {
    console.log('TODO not safe')
  }

  resolve(task.uid)
})

reducer.resolveOutput2 = (task, stream) => (dispatch) => new Promise((resolve, reject) => {
  const outerResolver = require('../resolvers.js')
  const chalk = require('chalk')

  stream.output = () => {
    console.log('  ' + chalk.italic('Resolving output'))

    // if (actionType === 'promise') {
    //   return stream._resolved
    // }

    // TODO this
    const customValidators = {
      fileResolvers: []
    }
    // if (props.validators && props.validators.file) {
    //   customValidators.fileResolvers = props.validators.file
    // }

    return outerResolver(task.output, task, 'end', customValidators)
  }

  // TODO actually dispatch
  resolve(task.uid)
})

function resolveInputHandler (state, action) {
  return Object.assign({}, state, {
    [action.uid]: Object.assign({}, state[action.uid], { resolvedInput: action.resolvedInput })
  })
}

reducer.resolveInput = (task) => (dispatch) => new Promise((resolve, reject) => {
  const outerResolver = require('../resolvers.js')

  // TODO handle from task lineage
  const resolvedInput = outerResolver(task.input, task, 'start')

  // let resolvedInput
  // if (passedOutput && !props.skipPassed) {
  //   // a. use values passed in from last Task OR
  //   resolvedInput = outerResolver(passedOutput, props, 'passed')
  //   console.log('Resolved from state:')
  //   console.log(resolvedInput)
  //
  //   // pop all DELETE_ME's
  //   // may break for more than one DELETE_ME..
  //   resolvedInput.forEach((item, i) => item === 'DELETE_ME' ? resolvedInput.splice(i, 1) : null)
  //
  //   console.log(resolvedInput)
  // } else {
  //   // b. resolve props to existing files
  //   console.log('resolving to fs')
  //   resolvedInput = outerResolver(input, props, 'start')
  // }

  dispatch({
    type: RESOLVE_INPUT,
    uid: task.uid,
    resolvedInput
  })

  resolve(task.uid)
})

function createActionHandler (state, action) {
  return Object.assign({}, state, {
    [action.uid]: Object.assign({}, state[action.uid], { action: action.action })
  })
}

reducer.createAction = (task, actionCreator) => (dispatch) => new Promise((resolve, reject) => {
  const resolvedProps = Object.assign({}, task, { input: task.resolvedInput })
  const action = actionCreator(resolvedProps)

  dispatch({
    type: CREATE_ACTION,
    uid: task.uid,
    action
  })

  resolve(task.uid)
})

reducer.setDuplex = (task, stream) => (dispatch) => new Promise((resolve, reject) => {
  const { isChildProcess } = require('../utils.js')
  const isStream = require('isstream')
  const { isWritable, isReadable, isDuplex } = isStream
  const action = task.action

  const actionType = (() => {
    if (isChildProcess(action)) {
      return 'process'
    } else if (action instanceof Promise) {
      return 'promise'
    } else if (typeof(action) === 'function') {
      // assume action returns function of signature fn(err, data)
      return 'curried-callback'
    } else if (isStream(action)) {
      return 'stream'
    } else {
      return 'unkown'
    }
  })()

  // 5. running the action
  switch(actionType) {
    case 'process':
      stream.setReadable(action)
      console.log(tab(1) + 'Creating a child process')
      break
    case 'promise':
      if (output && output.stream) {
      if (output.stream === 'object') {
        console.log(tab(1) + 'Creating a Readable for object from promise')

        // TODO do this earlier
        stream = duplexify.obj(null, null)
        stream.end()

        stream.setReadable(StreamFromPromise.obj(action))
        }
      } else {
        console.log(tab(1) + 'Creating a Readable for buffer/string from promise')
        stream.setReadable(StreamFromPromise(action))
      }
      break
    case 'curried-callback':
      stream.setReadable(new StreamFromCallback(action))
      break
    case 'stream':
      if (isReadable(action) && !isWritable(action)) {
        stream.setReadable(action)
        console.log('got a readable')
      } else if (isWritable(action) && !isReadable(action)) {
        stream.setWritable(action)
        console.log('got a writable')
      } else {
        // Duplex
        console.log('elsed to make a duplex')
        stream = action
      }
      // console.log(tab(1) + 'Creating a stream from stream')
      break
    default:
      // TODO this is a hack to make join work for now
      // stream = new Duplex()
      console.log('Bad action type')
  }


  // TODO dispatch something..

  resolve(task.uid)
})

reducer.catchTask = (stream) => (dispatch) => new Promise((resolve, reject) => {
  // IDEA: create stream._output on end, finish, events, then make
  // a close by calling stream.destroy()
  // (but wont work for child process)?
  stream.on('finish', function (chunk)  {
    console.log('  ' + 'Stream has finished')

    stream._output = stream.output()
    stream.destroy()

    if (chunk) {
      console.log('Last chunk ' + chunk)
    }
  })

  stream.on('close', () => {
    console.log('  ' + 'Child process finished')
  })

  stream.on('data', function(data) {
    if (actionType === 'promise') {
      this._resolved = data
    }
  })
})

module.exports = reducer
