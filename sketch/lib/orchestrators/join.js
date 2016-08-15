'use strict'

const _ = require('lodash')
const Promise = require('bluebird')

const { defaultCtx, mergeCtx } = require('../ctx')

function join(...tasks) {
  const accumulator = (currentCtx, task) => new Promise((resolve, reject) => {
    // Call next task with currentCtx
    task(_.noop, currentCtx).then((results) => {
      // Resolve to a new context with a ctx merge strategy
      const newCtx = mergeCtx('join')(currentCtx, results)
      resolve(newCtx)
    })
  })

  return (cb = _.noop, ctx = defaultCtx) => Promise.reduce(tasks, accumulator, ctx).asCallback(cb)
}

module.exports = join
