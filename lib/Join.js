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
      
      // TODO this will break when calling it on what join returns
      task(resolvedOutput)
        .on('task.done', (output) => { 
          if (i+1 !== tasks.length) {
            doTaskAndTryNext(i+1, tasks, output)
          } else if (i+1 === tasks.length) {
            emitter.emit('task.done', output)
          }
        })  
    }


    doTaskAndTryNext(0, tasks)

    // TODO return a valid Task
    // this kinda is, it follows the schema for task.done event
    return emitter
  }
}

module.exports = join
