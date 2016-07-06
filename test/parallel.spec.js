'use strict'

const { assert } = require('chai')

const task = require('../lib/Task.js')
const parallel = require('../lib/parallel.js')

describe('Parallel', function() {
  it('should run two tasks in parallel', function(done) {
    const task1 = task({ name: 'parallel 1' }, () => new Promise((resolve, reject) => {
      setTimeout(() => resolve('foo'), 1000)
    }))

    const task2 = task({ name: 'parallel 2' }, () => new Promise((resolve, reject) => {
      setTimeout(() => resolve('bar'), 500)
    }))

    parallel(task1, task2)()
      .on('task.done', (data) => {
        assert(data[0].toString() === 'bar')
        assert(data[1].toString() === 'foo') 
        
        done()
      })
  })
})
