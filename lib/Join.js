'use strict'

const { EventEmitter2 } = require('eventemitter2')
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
  return function() {
    const dup = duplexify(null, null)

    const doTaskAndTryNext = (i, tasks, resolvedOutput) => { 
      const task = tasks[i] 
     
      const finish = once((task) => {
        const output = task.output()

        if (i+1 !== tasks.length) {
          doTaskAndTryNext(i+1, tasks, output)
        } else if (i+1 === tasks.length) {
          // TODO change this
          dup.output = function() {
            return output
          }

          dup.destroy()
        }
      })


      // TODO this will break when calling it on what join returns
      task(resolvedOutput)
        // Assumes no chunk on close or finish
        // TODO handle a final chunk if it arrives
        // .on('data', (data) => console.log('Join received data: ', data))
        .on('end', function() {
          finish(this)
        })
        .on('close', function() {
          finish(this)
        })
        .on('finish', function() {
          finish(this)
        })
        .on('error', (err) => {
          console.log(err)
        })
    }

    doTaskAndTryNext(0, tasks)

    // TODO return a valid Task
    // has proper output, now needs input ability
    return dup
  }
}

module.exports = join
