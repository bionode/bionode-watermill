'use strict'

const { EventEmitter2 } = require('eventemitter2')

const { tab } = require('./utils.js')


/**
 * join
 *
 * @param ...tasks array of Tasks to join
 * @returns {undefined}
 */
function join(...tasks) {
  return function() {
    const emitter = new EventEmitter2()

    const doTaskAndTryNext = (i, tasks, resolvedOutput) => { 
      const task = tasks[i] 
     
      let finished = false
      const finish = (task) => {
        if (finished) {
          return
        }

        // Switch flag on
        finished = true

        const output = task.output()

        if (i+1 !== tasks.length) {
          doTaskAndTryNext(i+1, tasks, output)
        } else if (i+1 === tasks.length) {
          // TODO change this
          emitter.emit('task.done')
        }
      }


      // TODO this will break when calling it on what join returns
      task(resolvedOutput)
        // Assumes no chunk on close or finish
        // TODO handle a final chunk if it arrives
        // .on('end', function() {
        //   finish(this)
        // })
        .on('close', function() {
          finish(this)
        })
        .on('finish', function() {
          finish(this)
        })        
    }


    doTaskAndTryNext(0, tasks)

    // TODO return a valid Task
    // this kinda is, it follows the schema for task.done event
    return emitter
  }
}

module.exports = join
