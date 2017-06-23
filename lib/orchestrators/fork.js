'use strict'

const _ = require('lodash')
const Promise = require('bluebird')

const { defaultContext } = require('../constants/default-task-state.js')

function fork(...tasks) {
  const forkInvocator = (cb = _.noop, ctx = defaultContext) => {
    const taskToPromise = (task) => new Promise((resolve, reject) => {
      if (_.isArray(task.info)) {
        console.log('Fork created for: ', task.info.map(t => `${t.name} (${t.uid})`))
      }

      // Pass context into each forked task
      task(cb, ctx)
        .then(results => resolve(results))
        .catch(err => console.log('err: ', err))
    })

    return Promise.all(tasks.map(task => taskToPromise(task)))
      .then((results) => {
        return results
      }).asCallback(cb)
  }

  //forkInvocator.info = {
  //  uid: 'generic fork uid',
  //  name: 'generic fork name',
  // type: 'fork',
  //  tasks
  //}

  return forkInvocator
}

module.exports = fork
