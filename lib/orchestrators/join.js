'use strict'

const _ = require('lodash')
const Promise = require('bluebird')

const { defaultContext } = require('../constants/default-task-state.js')
const { mergeCtx } = require('../ctx')
const hash = require('../utils/hash.js')

function join(...tasks) {
  let uids = []
  const accumulator = (currentCtx, task) => new Promise((resolve, reject) => {
    // TODO better checking if task/join/junction/fork
    if (task.info) {
      console.log('Joining to task: ' + task.info.name)
      // Add this task to list of uids for this join
      uids.push(task.info.uid)
    }

    // Call next task with currentCtx
    task(_.noop, currentCtx).then((results) => {
      // Resolve to a new context with a ctx merge strategy
      const newCtx = mergeCtx('join')(currentCtx, results)
      resolve(newCtx)
    })
  })

  const joinInvocator = (cb = _.noop, ctx = defaultContext) =>
    Promise.reduce(tasks, accumulator, ctx)
      .then(results => Promise.resolve(Object.assign({}, {
        type: results.type,
        uid: hash(uids.join('')),
        tasks: uids,
        context: {
          trajectory: results.trajectory
        }
      })))
      .asCallback(cb)

  joinInvocator.info = tasks.map(task => task.info)

  return joinInvocator
}

module.exports = join
