'use strict'

const once = require('once')
const _ = require('lodash')
const { EventEmitter2 } = require('eventemitter2')
const Promise = require('bluebird')

const JOIN_FINISH = 'join.finish'

/**
 * join
 *
 * @param ...tasks array of Tasks to join
 * @returns {undefined}
 */
function join (...tasks) {
  return (trajectory = []) => {
    const emitter = new EventEmitter2()

    Promise.reduce(tasks, (currentTrajectory, task, i) => new Promise((resolve, reject) => {
      // Run the task and get its trajectory
      console.log('currentTrajectory: ', currentTrajectory)
      task(currentTrajectory).async
        .then((results) => {

          const { uid, params } = results
          let newTrajection = []

          // If it was a task (not a join or parallel)
          if (!_.isUndefined(uid) && !_.isUndefined(params)) {
            const paramsString = JSON.stringify(params)
            if (paramsString !== '{}') {
              newTrajection = [paramsString, uid]
            } else {
              newTrajection = [uid]
            }
          }

          newTrajection = results.trajectory.concat(newTrajection)

          console.log('newTrajection: ', newTrajection)

          resolve(newTrajection)
        })
    }), trajectory)
      .then((finalTrajectory) => {
        emitter.emit(JOIN_FINISH, {
          trajectory: finalTrajectory
        })
      })

    emitter.async = new Promise((resolve, reject) => {
      console.log('called join.async')
      emitter.on(JOIN_FINISH, (results) => resolve(results))
    })

    return emitter
  }
}

module.exports = join
