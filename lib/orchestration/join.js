'use strict'

const once = require('once')
const _ = require('lodash')
const { EventEmitter2 } = require('eventemitter2')

/**
 * join
 *
 * @param ...tasks array of Tasks to join
 * @returns {undefined}
 */
const join = (actions) => {
  return function join(...tasks) {
    return function(resolvedOutput) {
      const emitter = new EventEmitter2()

      const doTaskAndTryNext = (i, tasks, useCollection, trajection = [], uid) => {
        const task = tasks[i]

        const finish = once((task) => {
          const { uid, params } = task

          actions.addOutput(uid)
            .then((uid) => {
              let newTrajection
              const paramsString = JSON.stringify(params)
              if (paramsString !== '{}') {
                newTrajection = [paramsString, uid]
              } else {
                newTrajection = [uid]
              }

              if (trajection) {
                newTrajection = trajection.concat(newTrajection)
              }

              if (i+1 !== tasks.length) {
                doTaskAndTryNext(i+1, tasks, true, newTrajection, uid)
              } else if (i+1 === tasks.length) {
                emitter.emit('task.finish', {
                  resolvedOutput: 'collection'
                })
              }
            })
        })

        // console.log('trajection gonna be: ', trajection)

        task(useCollection, trajection)
          .on('task.finish', task => finish(task))
      }

      doTaskAndTryNext(0, tasks, resolvedOutput)

      // TODO return a valid Task
      // has proper output, now needs input ability
      return emitter
    }
  }
}

module.exports = join
