'use strict'

const _ = require('lodash')
const Promise = require('bluebird')

const { defaultContext } = require('../constants/default-task-state.js')
//const { mergeCtx } = require('../ctx')

const { addJunctionVertex } = require('../reducers/collection.js')
const hash = require('../utils/hash.js')
const orchestratorUidGenerator = require('./orchestrator-uid-generator')

const junction = (dispatch) => {
  return function (...tasks) {
    const junctionInvocator = (cb = _.noop, ctx = defaultContext) => {
      const taskToPromise = (task) => new Promise((resolve, reject) => {
        if (_.isArray(task.info.tasks)) {
          console.log('Junction created for: ', task.info)
        }

        // Pass context into each junctionized task
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
        dispatch(addJunctionVertex(uid, results))

        const junctionResults = {
          //type: 'junction',
          name: 'generic junction name',
          context: {
            trajectory: [uid]
          }
        }
        return junctionResults
      }).asCallback(cb)
    }

    junctionInvocator.info = {
      tasks: tasks,  // returns an array of tasks
      type: 'junction',  // used to check the type of orchestrator
      uid: hash(orchestratorUidGenerator(tasks).join(''))
    }

    //console.log('junctionInvocator')
    //console.log(junctionInvocator.info)

    return junctionInvocator
  }
}

module.exports = junction
