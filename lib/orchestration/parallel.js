'use strict'

const { EventEmitter2 } = require('eventemitter2')
const _ = require('lodash')

const PARALLEL_FINISH = 'parallel.finish'

function parallel (...tasks) {
  return (trajectory = []) => {
    const emitter = new EventEmitter2()
    const refTasks = new Array(tasks.length)

    const all = Promise.all(tasks.map((task, i) => new Promise((resolve, reject) => {
      const liveTask = task(trajectory)
      refTasks[i] = liveTask

      // liveTask may be task, join, or parallel
      // we can wildcard the event: *.finish
      // or use the async method
      liveTask.async.then((results) => {
        // TODO also emit the task/join/parallel.finish for internals here
        resolve(results.trajectory)
      })
    })))

    emitter.tasks = (cb) => {
      cb.call(null, refTasks)
      return emitter
    }

    emitter.async = new Promise((resolve, reject) => {
      emitter.on(PARALLEL_FINISH, (results) => resolve({
        trajectory: results
      }))
    })

    all.then((results) => {
      // Returns a flat concatenation of the trajectories from each task
      // Keeps order from original parallel call
      // e.g. parallel(A, join(B, C)) -> ['a', ['b', 'c']] -> ['a', 'b', 'c']
      // e.g. parallel(parallel(A, parallel(join(B, C), D), join(E, F))
      //      -> [ [ ['a'], [ ['b', 'c'], ['d'] ] ], ['e', 'f'] ]
      //      -> ['a', 'b', 'c', 'd', 'e', 'f']
      console.log('results: ', results)
      const flatResults = _.flattenDeep(results)
      // console.log(flatResults)
      emitter.emit(PARALLEL_FINISH, flatResults)
    })

    return emitter
  }
}

module.exports = parallel
