'use strict'

const stringify = require('json-stable-stringify')
const chalk = require('chalk')
const _ = require('lodash')
const endOfStream = require('end-of-stream')
const asyncDone = require('async-done')

const applicator = require('../utils/applicator.js')
const defaultConfig = require('./config.js').defaultState
const { hash, tab } = require('../utils')
const { defaultTask } = require('../constants/default-task-state.js')

// Actions
const CHECK_RESUMABLE = 'task/check-resumable'
const RESOLVE_OUTPUT = 'task/resolve-output'
const RESOLVE_INPUT = 'task/resolve-input'
const CREATE_ACTION = 'task/create-action'
const TASK_LOG = 'task/log'
const SET_VALIDATION = 'task/set-validation'
const ADD_TRAJECTION = 'task/add-trajection'

// Statuses
const {
  STATUS_CREATED,
  CHECKING_RESUMABLE,
  RESOLVING_INPUT
} = require('../constants/task-status-types.js')


/**
 * Initialize a task object using defaultTask and user passed in options.
 *
 * @param task {Object}
 * @returns a task object
 */
const initTask = (task) => {
  // null values NEED to be set
  // except for container, resolvedInput, resolvedOuput, operation
  // TODO a "deepObjectAssign", so that properties like params can be extended

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
    // TODO map action type to a function instead of tedius case statements
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
      break
    case SET_VALIDATION:
      return Object.assign({}, state, { validated: action.value })
      break
    case ADD_TRAJECTION:
      // const newTrajectory = state.trajectory.slice()
      // newTrajectory.push(action.trajection)
      return Object.assign({}, state, { trajectory: action.trajection })
      break
    default:
      // TODO second half of reconsider.
      // This is where we end up when proxied by tasks reducer
      // state will be the task from the tasks/create action ready to have a uid
      // created for it.
      return initTask(state)
  }
}


/**
 * Assigns task status to either RESOLVING_INPUT or CHECKING_RESUMABLE based on
 * `resume`.
 */
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


/**
 * Check if task has `resume` set to `on` of `off`
 *
 * TODO handle whether resolveInput or resolveOuput is called from here?
 */
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

  const { output, dir } = getState().tasks[uid]

  let opts = {}
  if (!_.isNull(dir)) {
    opts.cwd = dir
  }

  applicator(output, (item) => matchToFs(item, null, opts))
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
        .then(() => console.log('  ' + 'Wrote output(s) time+hash to watermill.json'))
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


/**
 * Simple assign `resolvedInput` to task state.
 *
 * TODO generic key assign handler
 */
function resolveInputHandler (state, { resolvedInput }) {
  return Object.assign({}, state, { resolvedInput })
}


/**
 * Resolve input.
 *
 * Takes `input` and runs it through a match-to-fs (TODO or otherwise)
 */
reducer.resolveInput = (uid, trajection) => (dispatch, getState) => new Promise((resolve, reject) => {
  const { input, dir } = getState().tasks[uid]

  if (_.isNull(input)) {
    dispatch({
      type: TASK_LOG,
      uid,
      content: tab(1) + 'Ignoring input as is null'
    })

    return resolve(uid)
  }

  dispatch({
    type: TASK_LOG,
    uid,
    content: tab(1) + 'Resolving input'
  })

  let opts = {}
  if (!_.isNull(dir)) {
    opts.cwd = dir
  }


  // TODO handle from task lineage
  const matchToFs = require('../matchers/match-to-fs.js')

  if (trajection.length === 0) {
    applicator(input, (item) => matchToFs(item, null, opts))
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
    const { trajectory } = getState().tasks[uid]
    console.log(`trajectory for ${uid}: `, trajectory)

    const minimatch = require('minimatch')
    const matchToCollection = (item) => new Promise((resolve, reject) => {
      let currentCollection = getState().collection

      const matchee = (path) => {
        if (minimatch(path.split('/').pop(), item)) {
          console.log(tab(2) + `${chalk.magenta(item)} matched to ${chalk.blue(path)}`)
          resolve(path)
        }
      }

      const traverse = (obj) => {
        if (_.isUndefined(obj)) {
          console.log('traverse got undefined, returning')
          return
        }

        if (_.isArray(obj)) {
          return obj.forEach(obj => traverse(obj))
        }

        if (_.isPlainObject(obj)) {
          return Object.keys(obj).forEach(key => traverse(obj[key]))
        }

        matchee(obj)
      }

      let didit = false
      if (trajectory.length > 0) {
        for (const node of trajectory.slice().reverse()) {
          traverse(currentCollection.vertexValue(node))
        }
      } else {
        console.log('Should only happen once')
        if (didit) {
          reject('only once')
        } else {
          didit = true
        }

        for (const [key, value] of currentCollection.vertices()) {
          // Assumes value !== undefined
          traverse(value)
        }
      }
    })


    applicator(input, matchToCollection)
      .then((resolvedInput) => {
        console.log('resolvedInput: ', resolvedInput)

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


/**
 * Simple assigns `operation` into task state.
 * TODO use generic key assign handler.
 */
function createActionHandler (state, { operation }) {
  return Object.assign({}, state, { operation })
}


/**
 * Create an `operation`. Resolve when it is done.
 *
 * Run operationCreator with resolved props (input resolved to fs) and set the
 * operation into the task state.
 */
reducer.createAction = (uid, operationCreator) => (dispatch, getState) => new Promise((resolve, reject) => {
  const task = getState().tasks[uid]
  // Replace `input` with `resolvedInput`
  const resolvedProps = Object.assign({}, task, { input: task.resolvedInput })
  // const operation = operationCreator(resolvedProps)

  // dispatch({
  //   type: CREATE_ACTION,
  //   uid,
  //   operation
  // })
  asyncDone(() => operationCreator(resolvedProps), (err, result) => {
    if (err) return reject(err)

    console.log('asyncDone result: ', result)
    resolve(uid)
  })
})

/**
 * Log.
 *
 * Simple logger action creator. The reducer handler will console.log the
 * content.
 *
 * TODO automatic indentation, custom logger
 */
reducer.taskLog = (uid, content) => (dispatch) => new Promise((resolve, reject) => {
  dispatch({
    type: TASK_LOG,
    uid,
    content
  })

  resolve(uid)
})

reducer.addTrajection = (uid, trajection) => (dispatch, getState) => new Promise((resolve, reject) => {
  dispatch({
    type: ADD_TRAJECTION,
    uid,
    trajection
  })

  resolve(uid)
})

module.exports = reducer
