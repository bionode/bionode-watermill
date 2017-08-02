'use strict'

const _ = require('lodash')
const Promise = require('bluebird')

const { defaultContext } = require('../constants/default-task-state.js')
const hash = require('../utils/hash.js')

const fork = () => {
  return function (...tasks) {
    const forkInvocator = (cb = _.noop, ctx = defaultContext) => {
      const taskToPromise = (task) => new Promise((resolve, reject) => {
        if (_.isArray(task.info.tasks)) {
          console.log('Fork created for: ', task.info)
        }

        // Pass context into each forked task
        task(cb, ctx)
          .then(results => resolve(results))
          .catch(err => console.log('err: ', err))
      })

      return Promise.all(tasks.map(task => taskToPromise(task)))
        .then((results) => {
          const uids = []
          for (const result of results) {
            if (result.uid) {
              uids.push(result.uid)
            } else {
              throw new Error('A uid was not found in the result')
            }
          }

          const uid = hash(uids.join(''))

          const forkResults = {
            name: 'generic fork name',
            context: {
              trajectory: [uid]
            }
          }
          return forkResults
        }).asCallback(cb)
    }

    forkInvocator.info = {
      tasks: tasks.map(task => task.info),  // returns an array of tasks
      type: 'fork'  // used to check the type of orchestrator
    }

    return forkInvocator
  }
}

module.exports = fork
