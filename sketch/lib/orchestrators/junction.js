'use strict'

const _ = require('lodash')
const Promise = require('bluebird')

const { defaultContext } = require('../constants/default-task-state.js')
const { mergeCtx } = require('../ctx')

const { addJunctionVertex } = require('../reducers/collection.js')

const junction = (dispatch) => {
  return function (...tasks) {
    const junctionInvocator = (cb = _.noop, ctx = defaultContext) => {
      const taskToPromise = (task) => new Promise((resolve, reject) => {
        if (_.isArray(task.info)) {
          console.log('Junction created for: ', task.info.map(t => `${t.name} (${t.uid})`))
        }

        // Pass context into each junctionized task
        task(cb, ctx)
        .then(results => resolve(results))
        .catch(err => console.log('err: ', err))
      })

      return Promise.all(tasks.map(task => taskToPromise(task)))
      .then((results) => {
        dispatch(addJunctionVertex(results))
        return Promise.resolve(results)
      })
      .then((results) => results.reduce(
        (sum, current) => Object.assign({}, sum, mergeCtx('junction')(sum, current)),
          defaultContext)
           )
           .then((results) => Object.assign({}, {
             type: 'junction',
             context: {
               trajectory: results.trajectory
             }
           }))
           .asCallback(cb)
    }

    return junctionInvocator
  }
}

module.exports = junction
