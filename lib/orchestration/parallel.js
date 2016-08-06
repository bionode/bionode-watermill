'use strict'

const duplexify = require('duplexify')
const { EventEmitter2 } = require('eventemitter2')
const once = require('once')

const noop = () => {}

function parallel(...tasks) {
  return function(trajection) {
    const emitter = new EventEmitter2()

    let total = tasks.length
    let completed = 0
    const results = []

    const tryFinish = () => {
      completed++

      if (completed === total) {
        emitter.emit('task.finish', {
          results
        })
      }
    }

    for (let task of tasks) {
      const finish = once((task) => {
        const output = task.resolvedOutput
        results.push(output)

        tryFinish()
      })

      task()
        //.on('data', noop) // force flowing mode
        .on('task.finish', (task) => finish(task) )
        .on('error', (err) => {
          console.log('ERROR: ', err)
        })
    }

    return emitter
  }
}

module.exports = parallel
