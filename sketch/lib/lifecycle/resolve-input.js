'use strict'

const _ = require('lodash')
const chalk = require('chalk')

const applicator = require('../utils/applicator.js')
const matchToFs = require('../matchers/match-to-fs.js')

/**
 * Lifecycle method: resolve input
 */

const dispatch = ({ content }) => console.log(content)
const tab = () => ''

const resolveInput = (taskState) => new Promise((resolve, reject) => {
  const { uid, input, dir, trajectory } = taskState

  if (_.isNull(input)) {
    console.log('Ignoring input as is null')

    return resolve(null)
  }

  console.log('resolving input')

  let opts = {}
  if (!_.isNull(dir)) {
    opts.cwd = dir
  }

  if (trajectory.length === 0) {
    applicator(input, (item) => matchToFs(item, null, opts))
      .then((resolvedInput) => {
        console.log(`Resolved ${chalk.magenta(JSON.stringify(resolvedInput, null, 2))} from ${chalk.blue(JSON.stringify(input, null, 2))}`)

        // dispatch({
        //   type: RESOLVE_INPUT,
        //   uid,
        //   resolvedInput
        // })
        resolve({ resolvedInput })
      })
      .catch(err => {
        console.log(err)
        reject('Input could not be resolved. Exiting early.')
      })
  } else {
    const { trajectory } = taskState
    console.log(`trajectory for ${uid}: `, trajectory)

    const minimatch = require('minimatch')
    const matchToCollection = (item) => new Promise((resolve, reject) => {
      let currentCollection = store.getState().collection

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

        // dispatch({
        //   type: RESOLVE_INPUT,
        //   uid,
        //   resolvedInput
        // })

        resolve({ resolvedInput })
      })
      .catch(err => console.log(err))
  }
})

module.exports = resolveInput
