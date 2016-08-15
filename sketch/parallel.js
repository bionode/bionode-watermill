'use strict'

const _ = require('lodash')
const Promise = require('bluebird')

const { defaultCtx, mergeCtx } = require('./ctx')

function parallel(...tasks) {
  return (cb = _.noop, ctx = defaultCtx) => Promise.all(tasks.map(task => new Promise((resolve, reject) => {
    // Pass context into each parallelized task
    task(_.noop, ctx).then(results => resolve(results))
  })))
  .then((results) => {
    return results.reduce((sum, current) => Object.assign({}, sum, mergeCtx('parallel')(sum, current)), defaultCtx)
  }).asCallback(cb)
}

module.exports = parallel
