'use strict'

const stringify = require('json-stable-stringify')
const defaultConfig = require('./config.js').defaultState
const { hash, tab } = require('../utils.js')
const chalk = require('chalk')
const _ = require('lodash')

// Actions
const CHECK_RESUMABLE = 'task/check-resumable'
const RESOLVE_OUTPUT = 'task/resolve-output'
const RESOLVE_INPUT = 'task/resolve-input'
const CREATE_ACTION = 'task/create-action'
const TASK_LOG = 'task/log'
const SET_VALIDATION = 'task/set-validation'

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
    validated: false,
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

// This reducer cannot have a default state is it is created based on the
// contents of the action. This is somewhat unconvential redux, but all actions
// beyond tasks/create act normal.
const reducer = (state, action) => {
  switch (action.type) {
    case CHECK_RESUMABLE:
      return checkResumableHandler(state, action)
      break
    case RESOLVE_OUTPUT:
      return resolveOutputHandler(state, action)
      break
    case RESOLVE_INPUT:
      return resolveInputHandler(state, action)
      break
    case CREATE_ACTION:
      return createActionHandler(state, action)
      break
    case TASK_LOG:
      // TODO store logs inside state
      // TODO automatic logging nesting
      console.log(action.content)
      return state
    case SET_VALIDATION:
      return Object.assign({}, state, { validated: action.value })
    default:
      // TODO second half of reconsider.
      // This is where we end up when proxied by tasks reducer
      // state will be the task from the tasks/create action ready to have a uid
      // created for it.
      return initTask(state)
  }
}

function checkResumableHandler(state, action) {
  const { resume } = state

  // TODO handle force-on and force-off as well
  if (resume === 'off') {
    // No resuming
    return Object.assign({}, state, { status: RESOLVING_INPUT })
  } else if (resume === 'on') {
    // Going to initiate resume attempt
    return Object.assign({}, state, { status: CHECKING_RESUMABLE })
  } else {
    throw new Error('task resume was an unexpected value or undefined')
    return state
  }
}

reducer.checkResumable = (uid) => (dispatch) => new Promise((resolve, reject) => {
  dispatch({
    type: CHECK_RESUMABLE,
    uid
  })

  resolve(uid)
})

function resolveOutputHandler(state, { resolvedOutput }) {
  return Object.assign({}, state, { resolvedOutput })
}

reducer.resolveOutput = (task) => (dispatch) => new Promise((resolve, reject) => {
  dispatch({
    type: TASK_LOG,
    uid: task.uid,
    content: tab(1) + 'Resolving output patterns to files'
  })

  const applicator = require('../applicator.js')
  const matchToFs = require('../matchers/match-to-fs.js')

  const { output } = task

  applicator(output, (item) => matchToFs(item, null))
  .then((resolvedOutput) => {
    dispatch({
      type: RESOLVE_OUTPUT,
      uid: task.uid,
      resolvedOutput
    })

    resolve(task.uid)
  })
  .catch(err => {
    // Not an actual error: one of the output files was not resolved
    console.log(tab(2) + err)
    return Promise.resolve()
  })
  .then(() => resolve(task.uid))
})

reducer.runValidators = (task, mode) => (dispatch) => new Promise((resolve, reject) => {
  dispatch({
    type: TASK_LOG,
    uid: task.uid,
    content: tab(1) + 'Running validations'
  })

  const { resolvedOutput } = task

  const applicator = require('../applicator.js')
  const {
    isNonNull,
    matchHash,
    matchTime
  } = require('../validators/')

  // TODO promise-ify other validators
  // validators.matchTime
  const validArr = [
    applicator(resolvedOutput, isNonNull)
  ]
  if (mode === 'before') {
    validArr.push(applicator(resolvedOutput, matchHash))
    validArr.push(applicator(resolvedOutput, matchTime))
  }

  const validations = Promise.all(validArr)

  // TODO dispatch something to store validations output

  validations
  .then(() => {
    if (task.status === CHECKING_RESUMABLE) {
      dispatch({
        type: TASK_LOG,
        uid: task.uid,
        content: tab(1) + 'All validations successful'
      })
    }

    dispatch({
      type: SET_VALIDATION,
      uid: task.uid,
      value: true
    })

    if (mode === 'after') {
      const checksummer = require('../checksummer.js')
      return applicator(resolvedOutput, checksummer)
        .then(() => console.log('  ' + 'Wrote output(s) time+hash to waterwheel.json'))
    }
  })
  .catch((err) => {
    dispatch({
      type: SET_VALIDATION,
      uid: task.uid,
      value: false
    })

    dispatch({
      type: TASK_LOG,
      uid: task.uid,
      content: tab(2) + err
    })

    dispatch({
      type: TASK_LOG,
      uid: task.uid,
      content: tab(1) + 'Could not resolve output before running'
    })
  })
  .then(() => resolve(task.uid))
})

function resolveInputHandler (state, { resolvedInput }) {
  return Object.assign({}, state, { resolvedInput })
}

reducer.resolveInput = (task, outputDump, trajectory) => (dispatch) => new Promise((resolve, reject) => {
  dispatch({
    type: TASK_LOG,
    uid: task.uid,
    content: tab(1) + 'Resolving input'
  })

  // TODO handle from task lineage
  const applicator = require('../applicator.js')
  const matchToFs = require('../matchers/match-to-fs.js')

  const { input } = task

  if (_.isNull(input)) {
    dispatch({
      type: TASK_LOG,
      uid: task.uid,
      content: tab(1) + 'Ignoring input as is null'
    })

    return resolve(task.uid)
  }

  if (!outputDump) {
    applicator(input, (item) => matchToFs(item, null))
      .then((resolvedInput) => {
        dispatch({
          type: TASK_LOG,
          uid: task.uid,
          content: tab(2) + `Resolved ${chalk.magenta(JSON.stringify(resolvedInput, null, 2))} from ${chalk.blue(JSON.stringify(input, null, 2))}`
        })

        dispatch({
          type: RESOLVE_INPUT,
          uid: task.uid,
          resolvedInput
        })

        resolve(task.uid)
      })
      .catch(err => {
        console.log(tab(2) + err)
        reject('Input could not be resolved. Exiting early.')
      })
  } else {
    const minimatch = require('minimatch')
    const matchToDump = (item) => new Promise((resolve, reject) => {
      let currentPlace = outputDump
      for (const trajection of trajectory) {
        currentPlace = currentPlace[trajection]

        for (const key in currentPlace.outputs) {
          const val = currentPlace.outputs[key]

          // TODO splitting to last piece b/c *.foo patterns do not work
          // change this, or figure out edge cases
          const doer = (val) => new Promise((resolve, reject) => {
            console.log(val, item)
            if (minimatch(val.split('/').pop(), item)) {
              console.log(tab(2) + `${chalk.magenta(item)} matched to ${chalk.blue(val)}`)
              resolve(val)
            } else {
              // this feels like it will happen too much...
              reject()
            }
          })

          applicator(val, doer)
            .then(val => resolve(val))
            .catch(err => reject('Did not make a match: ', err))
        }
      }
    })

    applicator(input, matchToDump)
      .then((resolvedInput) => {
        dispatch({
          type: RESOLVE_INPUT,
          uid: task.uid,
          resolvedInput
        })


        resolve(task.uid)
      })
      .catch(err => console.log(err))
  }
})

// TODO better name than 'action'
function createActionHandler (state, { action }) {
  return Object.assign({}, state, { action })
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
        // console.log('got a readable')
      } else if (isWritable(action) && !isReadable(action)) {
        stream.setWritable(action)
        // console.log('got a writable')
      } else {
        // Duplex
        console.log('elsed to make a duplex')

        // TODO sketchy
        // const duplexify = require('duplexify')
        // stream = duplexify.obj(null, null)
        // console.log(action)
        stream.setReadable(action)
        // Initiale flowing mode
        action.on('data', () => {})
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

reducer.catchTask = (task, stream) => (dispatch) => new Promise((resolve, reject) => {
  // IDEA: create stream._output on end, finish, events, then make
  // a close by calling stream.destroy()
  // (but wont work for child process)?
  stream.on('finish', function (chunk)  {
    console.log('  ' + 'Stream has finished')

    // stream._output = task.resolvedOutput
    // stream.destroy()
    resolve(task.uid)

    if (chunk) {
      console.log('Last chunk ' + chunk)
    }
  })

  stream.on('end', function(chunk) {
    console.log('  ' + 'Received end internally')
    resolve(task.uid)
  })

  stream.on('close', () => {
    console.log('  '+ 'Received close internally')
    resolve(task.uid)
  })

  stream.on('data', function(data) {
    if (actionType === 'promise') {
      this._resolved = data
    }
  })

  // TODO resolve(uid)
})

reducer.taskLog = (uid, content) => (dispatch) => new Promise((resolve, reject) => {
  dispatch({
    type: TASK_LOG,
    uid,
    content
  })

  resolve(uid)
})

module.exports = reducer
