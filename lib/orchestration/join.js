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
          const { uid } = task

          actions.addOutput(uid)
            .then((uid) => {
              let newTrajection
              if (trajection) {
                newTrajection = trajection.slice()
                newTrajection.push(uid)
              } else {
                newTrajection = [uid]
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
