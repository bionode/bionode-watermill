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
    return function(trajection) {
      const emitter = new EventEmitter2()
      const ids = new Array(tasks.length)

      const doTaskAndTryNext = (i, tasks, trajection = [], uid) => {
        const task = tasks[i]

        const finish = once((task) => {
          console.log('trajection from just inside finish: ', task.trajection)
          const { uid, params } = task

          ids[i] = uid

          if (!_.isUndefined(uid)) {
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

              console.log('newTrajection: ', newTrajection)
              if (i+1 !== tasks.length) {
                doTaskAndTryNext(i+1, tasks, newTrajection, uid)
              } else if (i+1 === tasks.length) {
                console.log('emitting a finish from: ', uid)
                emitter.emit('task.finish', {
                  trajection: newTrajection,
                  name: 'joiner finish',
                  joined: ids
                })
              }
            })
          } else {
            console.log('undefined uid, must be from a join join')
            console.log('task.trajection: ', task.trajection)
            emitter.emit('task.finish', {
              trajection: task.trajection,
              name: 'end of joiner'
            })
          }
        })


        console.log('calling a task, with trajection: ', trajection)
        task(trajection)
          .on('task.finish', task => finish(task))
      }

      console.log('tasks: ', tasks)
      doTaskAndTryNext(0, tasks, trajection)

      // TODO return a valid Task
      // has proper output, now needs input ability
      return emitter
    }
  }
}

module.exports = join
