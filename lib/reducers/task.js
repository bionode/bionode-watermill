'use strict'

const stringify = require('json-stable-stringify')
const chalk = require('chalk')
const _ = require('lodash')
const endOfStream = require('end-of-stream')

const applicator = require('../utils/applicator.js')
const defaultConfig = require('./config.js').defaultState
const { hash, tab } = require('../utils')

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
    params: {},
    trajectory: []
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


/**
 * Resolve output.
 *
 * Takes `output` from task object, returns an object with same structure but
 * with `matchToFs` (OR TODO another matcher passed in) ran on each leaf of the
 * object.
 *
 * Dispatches resolve output on success, otherwise only resolves to task uid.
 * It will be another action creator's duty to check if `resolvedOutput` is
 * present and act accordingly.
 */
reducer.resolveOutput = (uid) => (dispatch, getState) => new Promise((resolve, reject) => {
  dispatch({
    type: TASK_LOG,
    uid,
    content: tab(1) + 'Resolving output patterns to files'
  })

  const matchToFs = require('../matchers/match-to-fs.js')

  const { output } = getState().tasks[uid]

  applicator(output, (item) => matchToFs(item, null))
  .then((resolvedOutput) => {
    dispatch({
      type: RESOLVE_OUTPUT,
      uid,
      resolvedOutput
    })

    resolve(uid)
  })
  .catch(err => {
    // Not an actual error: one of the output files was not resolved
    dispatch({
      type: TASK_LOG,
      uid,
      content: tab(2) + err
    })

    return Promise.resolve()
  })
  .then(() => resolve(uid))
})


/**
 * Run validations.
 *
 * @param mode {String} `before` or `after` TODO comes from task status
 *
 * Takes `resolvedOuput` from task and runs validators over each leaf.
 */
reducer.runValidators = (uid, mode) => (dispatch, getState) => new Promise((resolve, reject) => {
  dispatch({
    type: TASK_LOG,
    uid,
    content: tab(1) + 'Running validations'
  })

  const { resolvedOutput, status } = getState().tasks[uid]

  const {
    isNonNull,
    matchHash,
    matchTime
  } = require('../validators')

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
    if (status === CHECKING_RESUMABLE) {
      dispatch({
        type: TASK_LOG,
        uid,
        content: tab(1) + 'All validations successful'
      })
    }

    dispatch({
      type: SET_VALIDATION,
      uid,
      value: true
    })

    if (mode === 'after') {
      const checksummer = require('../post-validators/checksummer.js')
      return applicator(resolvedOutput, checksummer)
        .then(() => console.log('  ' + 'Wrote output(s) time+hash to waterwheel.json'))
    }
  })
  .catch((err) => {
    dispatch({
      type: SET_VALIDATION,
      uid,
      value: false
    })

    dispatch({
      type: TASK_LOG,
      uid,
      content: tab(2) + err
    })

    dispatch({
      type: TASK_LOG,
      uid,
      content: tab(1) + 'Could not resolve output before running'
    })
  })
  .then(() => resolve(uid))
})

function resolveInputHandler (state, { resolvedInput }) {
  return Object.assign({}, state, { resolvedInput })
}

reducer.resolveInput = (uid, outputDump, trajectory) => (dispatch, getState) => new Promise((resolve, reject) => {
  dispatch({
    type: TASK_LOG,
    uid,
    content: tab(1) + 'Resolving input'
  })

  // TODO handle from task lineage
  const matchToFs = require('../matchers/match-to-fs.js')

  const { input } = getState().tasks[uid]

  if (_.isNull(input)) {
    dispatch({
      type: TASK_LOG,
      uid,
      content: tab(1) + 'Ignoring input as is null'
    })

    return resolve(uid)
  }

  if (!outputDump) {
    applicator(input, (item) => matchToFs(item, null))
      .then((resolvedInput) => {
        dispatch({
          type: TASK_LOG,
          uid,
          content: tab(2) + `Resolved ${chalk.magenta(JSON.stringify(resolvedInput, null, 2))} from ${chalk.blue(JSON.stringify(input, null, 2))}`
        })

        dispatch({
          type: RESOLVE_INPUT,
          uid,
          resolvedInput
        })

        resolve(uid)
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
              reject(`${item} did not match ${val}`)
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
          uid,
          resolvedInput
        })


        resolve(uid)
      })
      .catch(err => console.log(err))
  }
})

// TODO better name than 'action'
function createActionHandler (state, { action }) {
  return Object.assign({}, state, { action })
}

reducer.createAction = (uid, operationCreator) => (dispatch, getState) => new Promise((resolve, reject) => {
  const task = getState().tasks[uid]
  // Replace `input` with `resolvedInput`
  const resolvedProps = Object.assign({}, task, { input: task.resolvedInput })
  const operation = operationCreator(resolvedProps)

  dispatch({
    type: CREATE_ACTION,
    uid,
    action: operation
  })

  resolve(uid)
})

reducer.setDuplex = (uid, stream) => (dispatch, getState) => new Promise((resolve, reject) => {
  const { isChildProcess } = require('../utils')
  const isStream = require('isstream')
  const { isWritable, isReadable, isDuplex } = isStream

  const { action } = getState().tasks[uid]

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

  // readable, writable, or duplex so that end-of-stream can work
  let streamType

  // 5. running the action
  switch(actionType) {
    case 'process':
      stream.setReadable(action)
      // TODO dispatch streamType into store
      streamType = 'readable'
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
      streamType = 'readable'
      break
    case 'curried-callback':
      stream.setReadable(new StreamFromCallback(action))
      streamType = 'readable'
      break
    case 'stream':
      if (isReadable(action) && !isWritable(action)) {
        stream.setReadable(action)
        streamType = 'readable'
        console.log('got a readable')
      } else if (isWritable(action) && !isReadable(action)) {
        stream.setWritable(action)
        streamType = 'writable'
        console.log('got a writable')
      } else {
        // Duplex
        console.log('elsed to make a duplex')

        // TODO sketchy
        // const duplexify = require('duplexify')
        // stream = duplexify.obj(null, null)
        // console.log(action)
        stream.setReadable(action)
        stream.setWritable(action)
        streamType = 'duplex'
        // Initiale flowing/legacy mode
        // TODO not this
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
  // TODO resolve only with uid
  resolve([uid, streamType])
})

reducer.catchStream = (uid, stream, streamType) => (dispatch, getState) => new Promise((resolve, reject) => {
  let opts = {}
  if (streamType === 'readable') {
    opts.writable = false
  } else if (streamType === 'writable') {
    opts.readable = false
  }

  // TODO how to handle watching readable and/or writable portions of duplex?
  endOfStream(stream, opts, (err) => {
    if (err) return reject(err)

    resolve(uid)
  })
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
