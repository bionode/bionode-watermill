const { assert } = require('chai')

let task, shell

describe('Task', function() {
  beforeEach(function() {
    const waterwheelInstance = require('../')()
    task = waterwheelInstance.task
    join = waterwheelInstance.join
    shell = waterwheelInstance.shell
  })

  describe('on finish event', function() {
    it('should have an empty array as trajectory if none is passed', function(done) {
      const myTask = task({
        input: null,
        output: 'foo.txt'
      }, () => shell('echo foo > foo.txt'))

      myTask()
        .on('task.finish', (results) => {
          const { trajectory } = results

          assert.deepEqual(trajectory, [])
          done()
        })
    })

    it.skip('should have the passed in trajectory', function(done) {
      const myTask = task({
        input: null,
        output: 'foo.txt'
      }, () => shell('echo foo > foo.txt'))

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
      const myTask = task({
        input: null,
        output: 'foo.txt'
      }, () => shell('echo foo > foo.txt'))().async.then((results) => {
        const { trajectory } = results

        assert.deepEqual(trajectory, [])
        done()
      })
    })

    it.skip('should have the passed in trajectory', function(done) {
      task()(['a', 'b']).async.then((results) => {
        const { trajectory } = results

        assert.deepEqual(trajectory, ['a', 'b'])
        done()
      })
    })
  })
})

describe('Join', function() {
  it.skip('should build a trajectory', function(done) {
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
  it.skip('should run all tasks', function(done) {
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
