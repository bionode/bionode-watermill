/*
 * This is a simple implementation of the task orchestration utilities provided
 * by bionode-waterwheel. The implementation is complemented by a test suite.
 * A fully realized implementation shall satisfy the test suite, regardless of
 * internal details. This means alternative implementations can be used as long
 * as they implement the abstract API documented here and pass this test suite.
 */

const { EventEmitter2 } = require('eventemitter2')
const Promise = require('bluebird')

const TASK_FINISH = 'task.finish'
const JOIN_FINISH = 'join.finish'
const PARALLEL_FINISH = 'parallel.finish'

let alphabet = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n', 'o', 'p', 'q', 'r', 's', 't', 'u', 'v', 'w', 'x', 'y', 'z']

/**
 * === TASK ===
 * A **task** is the fundamental **unit** for building pipelines.
 *
 * @params props {Object} describes I/O, settings, etc.
 * @returns {anonymous function} the "task starter"
 *
 * @params trajectory {Array} describes location in pipeline
 * @returns {EventEmitter} emits `task.finish` with payload including trajectory
 */
const task = (props) => (trajectory = []) => {
  // NOTE props is not required for a minimal implementation
  const emitter = new EventEmitter2()

  const results = {
    trajectory
  }

  emitter.async = new Promise((resolve, reject) => {
    emitter.on(TASK_FINISH, resolve)
  })

  setTimeout(() => emitter.emit(TASK_FINISH, results))

  return emitter
}

/**
 * === JOIN ===
 *
 * task1 -> task2
 */
function join (...tasks) {
  return (trajectory = []) => {
    const emitter = new EventEmitter2()

    emitter.async = new Promise((resolve, reject) => {
      emitter.on(JOIN_FINISH, resolve)
    })

    Promise.reduce(tasks, function(currentTrajectory, task, i) {
      return new Promise((resolve, reject) => {
        // Run the task and get its trajectory
        task([alphabet[i]]).async.then((results) => {
          resolve(currentTrajectory.concat(results.trajectory))
        })
      })
    }, trajectory).then((finalTrajectory) => {
      emitter.emit('join.finish', {
        trajectory: finalTrajectory
      })
    })

    return emitter
  }
}

function parallel (...tasks) {
  return (trajectory = []) => {
    const emitter = new EventEmitter2()
    const refTasks = new Array(tasks.length)

    const all = Promise.all(tasks.map((task, i) => new Promise((resolve, reject) => {
      const liveTask = task(trajectory)
      refTasks[i] = liveTask

      liveTask.on('task.finish', (results) => {
        // emitter.emit('task.finish', results)
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
      emitter.emit(PARALLEL_FINISH, results)
    })

    return emitter
  }
}

// === TEST SUITE ===
const { assert } = require('chai')

describe('Task', function() {
  describe('on finish event', function() {
    it('should have an empty array as trajectory if none is passed', function(done) {
      const myTask = task()

      myTask()
        .on('task.finish', (results) => {
          const { trajectory } = results

          assert.deepEqual(trajectory, [])
          done()
        })
    })

    it('should have the passed in trajectory', function(done) {
      const myTask = task()

      myTask(['a', 'b'])
        .on('task.finish', (results) => {
          const { trajectory } = results

          assert.deepEqual(trajectory, ['a', 'b'])
          done()
        })
    })
  })

  describe('promise method', function() {
    it('should have an empty array as trajectory if none is passed', function(done) {
      task()().async.then((results) => {
        const { trajectory } = results

        assert.deepEqual(trajectory, [])
        done()
      })
    })

    it('should have the passed in trajectory', function(done) {
      task()(['a', 'b']).async.then((results) => {
        const { trajectory } = results

        assert.deepEqual(trajectory, ['a', 'b'])
        done()
      })
    })
  })
})

describe('Join', function() {
  it('should build a trajectory', function(done) {
    const taskA = task()
    const taskB = task()
    const taskC = task()

    join(taskA, taskB, taskC)().async.then((results) => {
      const { trajectory } = results
      assert.deepEqual(trajectory, ['a', 'b', 'c'])
      done()
    })
  })
})

describe('Parallel', function() {
  it('should run all tasks', function(done) {
    const taskA = task()
    const taskB = task()

    parallel(taskA, taskB)()
      .async.then((results) => {
        const { trajectory } = results
        assert.deepEqual(trajectory, [[], []])
        done()
      })
  })
})
