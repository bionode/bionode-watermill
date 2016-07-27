'use strict'

const once = require('once')
const duplexify = require('duplexify')

const { tab } = require('./utils.js')

/**
 * join
 *
 * @param ...tasks array of Tasks to join
 * @returns {undefined}
 */
function join(...tasks) {
  let store = []

  return function(resolvedOutput) {
    const dup = duplexify(null, null)

    const doTaskAndTryNext = (i, tasks, resolvedOutput) => {
      const task = tasks[i]

      const finish = once((task) => {
        const output = task._output

        // TODO this better
        if (Array.isArray(output)) {
          // will fail if it is an array of arrays... just a hotfix
          output.forEach(item => store.push(item))
        } else {
          store.push(output)
        }

        if (i+1 !== tasks.length) {
          // console.log('store:')
          // console.log(store)
          doTaskAndTryNext(i+1, tasks, store)
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
      task(resolvedOutput)
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
          finish(this)
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
