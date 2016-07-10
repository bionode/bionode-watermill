'use strict'

const duplexify = require('duplexify')
const once = require('once')

const noop = () => {}

function parallel(...tasks) {
  return function() {
    const dup = duplexify(null, null)

    let total = tasks.length
    let completed = 0
    const results = []

    const tryFinish = () => {
      completed++

      if (completed === total) {
        dup.output = () => {
          return results
        }
        dup.destroy()
      }
    }

    for (let task of tasks) {
      const finish = once((task) => {
        const output = task.output()
        results.push(output)

        tryFinish()
      })

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

    return dup
  }
}

module.exports = parallel
