'use strict'

const _ = require('lodash')
const Promise = require('bluebird')

const { defaultContext } = require('../constants/default-task-state.js')
const hash = require('../utils/hash.js')
const orchestratorUidGenerator = require('./orchestrator-uid-generator')

const fork = (dispatch) => {
  return function (...tasks) {
    const forkInvocator = (cb = _.noop, ctx = defaultContext) => {
      const taskToPromise = (task) => new Promise((resolve, reject) => {
        if (_.isArray(task.info.tasks)) {
          console.log('Fork created for: ', task.info)
        }
        console.log('fork ctx: ', ctx)
         // Pass context into each forked task
         task(cb, ctx)
           .then(results => resolve(results))
           .catch(err => console.log('err: ', err))
      })

      return Promise.all(tasks.map(task => taskToPromise(task)))
        .then((results) => {
        console.log("fork results: ", results)
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
          console.log("return fork: ", forkResults)
          return forkResults
        }).asCallback(cb)
    }

    forkInvocator.info = {
      tasks: tasks,  // returns an array of tasks
      type: 'fork',  // used to check the type of orchestrator
      uid: hash(orchestratorUidGenerator(tasks).join(''))
    }

    //console.log("forkInvocator")
    //console.log(forkInvocator.info)

    return forkInvocator
  }

}

module.exports = fork
