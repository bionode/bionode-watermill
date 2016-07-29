'use strict'

const once = require('once')
const duplexify = require('duplexify')
const _ = require('lodash')

const { tab } = require('../utils')

/**
 * join
 *
 * @param ...tasks array of Tasks to join
 * @returns {undefined}
 */
function join(...tasks) {
  let store = {}

  return function(resolvedOutput) {
    const dup = duplexify(null, null)

    const doTaskAndTryNext = (i, tasks, resolvedOutput, trajectory) => {
      const task = tasks[i]

      const finish = once((task) => {
        const output = task._output

        // TODO this better

        if (trajectory) {
          task.stateWormhole.trajectory = trajectory
        }

        const paramsStr = JSON.stringify(task.stateWormhole.params)
        task.stateWormhole.trajectory.push(paramsStr)
        console.log('trajectory: ', task.stateWormhole.trajectory)

        if (Object.keys(store).length === 0) {
          store = { [paramsStr]: {} }
        }

        let currentPlace = store
        // for (const trajection of task.stateWormhole.trajectory) {
        task.stateWormhole.trajectory.forEach( (trajection, i) => {
          if (currentPlace[trajection]) {
            currentPlace = currentPlace[trajection]
          } else {
            currentPlace[trajection] = {}
          }

          if (i === task.stateWormhole.trajectory.length - 1) {
            if (paramsStr !== trajection) {
              // Same params string as last trajection
              currentPlace[paramsStr] = {}
              currentPlace = currentPlace[paramsStr]
            }

            console.log(currentPlace)

            Object.assign(currentPlace, { [task.stateWormhole.uid]:output })
          }

        })

        console.log('store:', JSON.stringify(store, null, 2))

        if (i+1 !== tasks.length) {
          // console.log('store:')
          // console.log(store)
          doTaskAndTryNext(i+1, tasks, store, task.stateWormhole.trajectory)
        } else if (i+1 === tasks.length) {
          // TODO change this
          dup.output = function() {
            return store
          }
          dup._output = store

          dup.destroy()
          dup.emit('destroy')
        }
      })


      // TODO this will break when calling it on what join returns
      task(resolvedOutput, trajectory)
        // Assumes no chunk on close or finish
        // TODO handle a final chunk if it arrives
        // .on('data', (data) => console.log('Join received data: ', data))
        // .on('end', function() {
        //   // console.log('Join received: end')
        //   finish(this)
        // })
        // .on('close', function() {
        //   // console.log('Join received: close')
        //   finish(this)
        // })
        // .on('finish', function() {
        //   // console.log('Join received: finish')
        //   finish(this)
        // })
        .on('destroy', function() {
          finish(this, task)
        })
        .on('error', (err) => {
          console.log(err)
        })
    }

    doTaskAndTryNext(0, tasks, resolvedOutput)

    // TODO return a valid Task
    // has proper output, now needs input ability
    return dup
  }
}

module.exports = join
