'use strict'

const _ = require('lodash')
const chalk = require('chalk')
const fs = require('fs')

const applicator = require('../utils/applicator.js')
const matchToFs = require('../matchers/match-to-fs.js')

/**
 * Lifecycle method: resolve input
 */

//const dispatch = ({ content }) => console.log(content)
const tab = () => ''

const resolveInput = (taskState, DAG, logger) => new Promise((resolve, reject) => {
  const { uid, input, dir } = taskState
  const miniUid = uid.substring(0, 7)
  const { trajectory } = taskState.context
  logger.emit('log', `Resolving input for ${uid.substring(0,7)}`)

  if (_.isNull(input)) {
    logger.emit('log', 'Ignoring input as is null', 1)

    logger.emit('log', `Successfully resolved input for ${miniUid}`)
    return resolve(null)
  }

  let opts = {}
  // TODO use dir
  // if (!_.isNull(dir)) {
  //   opts.cwd = dir
  // }

  if (trajectory.length === 0) {
    applicator(input, (item) => matchToFs(item, null, opts))
      .then((resolvedInput) => {
        console.log(`Resolved ${chalk.magenta(JSON.stringify(resolvedInput, null, 2))} from ${chalk.blue(JSON.stringify(input, null, 2))}`)

        logger.emit('log', `Successfully resolved input for ${uid}`)

        resolve({ resolvedInput })
      })
      .catch(err => {
        console.log(err)
        reject('Input could not be resolved. Exiting early.')
      })
  } else {
    console.log(`trajectory for ${uid}: `, trajectory)

    const minimatch = require('minimatch')
    const matchToCollection = (item) => new Promise((resolve, reject) => {
      let currentCollection = DAG
      console.log('currentCollection: ', JSON.stringify(require('../reducers/collection.js').jsonifyGraph(currentCollection), null, 2))

      const matchee = (path) => {
        if (minimatch(path.split('/').pop(), item)) {
          console.log(tab(2) + `${chalk.magenta(item)} matched to ${chalk.blue(path)}`)
          return resolve(path)
        } else if (minimatch(path, item)) {
          console.log(tab(2) + `${chalk.magenta(item)} matched to ${chalk.blue(path)}`)
          return resolve(path)
        } else {
          // this checks wether the input is on current working directory
          // and tries to match it to input pattern
          fs.readdir(process.cwd(), (err, files) => {
            if (err) {
              throw new Error(err)
            }
            files.forEach(file => {
              if (minimatch(file, item)) {
                console.log(tab(2) + `${chalk.magenta(item)} matched to ${chalk.blue(file)}`)
                return resolve(file)
              }
            })
          })
        }
      }

      const traverse = (obj) => {
        if (_.isUndefined(obj)) {
          logger.emit('log', 'traverse got undefined, returning')
          return
        }

        if (_.isArray(obj)) {
          if (obj.length === 0) {
            return reject('No values to match to at DAG node')
          }
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

        resolve({ resolvedInput })
      })
      // .catch(err => console.log('err: ', err))
      .catch(err => reject(err))
  }
})

module.exports = resolveInput
