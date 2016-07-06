'use strict'

const { EventEmitter2 } = require('eventemitter2')
const noop = () => {}

function parallel(...tasks) {
  return function() {
    const emitter = new EventEmitter2()

    let total = tasks.length
    let completed = 0
    const results = []

    const tryFinish = () => {
      completed++

      if (completed === total) {
        emitter.emit('task.done', results)
      }
    }


    for (let task of tasks) {
      let finished = false
      const finish = (task) => {
        if (finished) {
          return
        }

        // Switch flag on
        finished = true

        const output = task.output()
        results.push(output)

        tryFinish()
      }


      task()
        .on('data', noop) // force flowing mode
        .on('end', function() {
          finish(this)
        })
        .on('finish', function() {
          finish(this)
        })
        .on('error', (err) => {
          console.log('ERROR: ', err)
        })
    }

    return emitter
  }
}

module.exports = parallel
