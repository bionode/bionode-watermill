'use strict'

const fs = require('fs')
const path = require('path')
const { assert } = require('chai')

let task, parallel

const twd = path.resolve(__dirname, 'files', 'parallel')

const delayedEnd = (str, time) => task({
  name: `Return ${str} after ${time}ms`
}, () => new Promise((resolve, reject) => {
  setTimeout(() => resolve(str), time)
}))

const delayedWrite = (str, time) => task({
  input: null,
  output: `parallel.spec.${str}.txt`,
  params: { data: str },
  dir: twd,
  name: `Readable to end after ${time}ms`
}, ({ params, output, dir }) => {
  const rs = require('readable-stream').Readable({ read: () => {} })
  setTimeout(() => {
    rs.push(params.data)
    rs.push(null)
  }, time)

  return rs.pipe(fs.createWriteStream(dir + '/' + output))
})

describe('Parallel', function() {
  beforeEach(function() {
    const waterwheelInstance = require('../')()
    task = waterwheelInstance.task
    parallel = waterwheelInstance.parallel
  })

  it('should run two tasks in parallel', function(done) {
    const task1 = delayedWrite('foo', 1000)
    const task2 = delayedWrite('bar', 500)

    parallel(task1, task2)().async.then((results) => {
      console.log('results: ', results)
      done()
    })
  })

  it.skip('should parallel parallel', function(done) {
    const taskA1 = delayedWrite('A1', 100)
    const taskA2 = delayedWrite('A2', 500)
    const taskB1 = delayedWrite('B1', 200)
    const taskB2 = delayedWrite('B2', 400)

    const p1 = parallel(taskA1, taskA2)
    const p2 = parallel(taskB1, taskB2)

    parallel(p1, p2)().async.then((results) => {
      const data = results.trajectory

      assert(data[0][0].toString() === 'B1')
      assert(data[0][1].toString() === 'B2')
      assert(data[1][0].toString() === 'A1')
      assert(data[1][1].toString() === 'A2')

      done()
    })
  })

  it.skip('should arrive in proper order', function(done) {
    const taskA1 = delayedResolve('A1', 100)
    const taskA2 = delayedResolve('A2', 500)
    const taskB1 = delayedResolve('B1', 200)
    const taskB2 = delayedResolve('B2', 400)

    parallel(taskA1, taskA2, taskB1, taskB2)()
      .on('close', function() {
        const data = this.output()

        assert(data[0].toString() === 'A1')
        assert(data[1].toString() === 'B1')
        assert(data[2].toString() === 'B2')
        assert(data[3].toString() === 'A2')

        done()
      })
  })
})
