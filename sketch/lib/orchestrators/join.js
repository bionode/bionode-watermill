'use strict'

const _ = require('lodash')
const Promise = require('bluebird')

const { defaultContext } = require('../constants/default-task-state.js')
const { mergeCtx } = require('../ctx')

function join(...tasks) {
  const accumulator = (currentCtx, task) => new Promise((resolve, reject) => {
    console.log('Joining to task: ' + task.info.name)

    // Call next task with currentCtx
    task(_.noop, currentCtx).then((results) => {
      // Resolve to a new context with a ctx merge strategy
      const newCtx = mergeCtx('join')(currentCtx, results)
      resolve(newCtx)
    })
  })

  return (cb = _.noop, ctx = defaultContext) =>
    Promise.reduce(tasks, accumulator, ctx)
      .then(results => Promise.resolve(Object.assign({}, {
        type: results.type,
        context: {
          trajectory: results.trajectory
        }
      })))
      .asCallback(cb)
}

module.exports = join
