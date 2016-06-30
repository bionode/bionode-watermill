'use strict'

const { tab } = require('./utils.js')

function join(...tasks) {
  const doTaskAndTryNext = (i, tasks) => { 
    const task = tasks[i]

    task()
      .on('task.done', (output) => {
        console.log('Task finished with output:')
        console.log(tab(1) + output)
  
        if (i+1 !== tasks.length) {
          doTaskAndTryNext(i+1, tasks)
        }
      })  
  }


  doTaskAndTryNext(0, tasks)

  // for (let i=0; i < tasks.length - 1; i++) {
  //   const task = tasks[i]
  //   
  //   task()
  // } 
}

module.exports = join
