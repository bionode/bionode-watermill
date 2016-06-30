'use strict'

const { tab } = require('./utils.js')


/**
 * join
 *
 * @param ...tasks array of Tasks to join
 * @returns {undefined}
 */
function join(...tasks) {
  const doTaskAndTryNext = (i, tasks, resolvedOutput) => { 
    const task = tasks[i] 

    task(resolvedOutput)
      .on('task.done', (output) => { 
        if (i+1 !== tasks.length) {
          doTaskAndTryNext(i+1, tasks, output)
        }
      })  
  }


  doTaskAndTryNext(0, tasks)

  // TODO return a valid Task
}

module.exports = join
