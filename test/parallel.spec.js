'use strict'

const fs = require('fs')
const path = require('path')
const { assert } = require('chai')

const { task, parallel } = require('../')

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

//delayedWrite('foo', 1000)()
//.on('task.finish', (task) => {
//  console.log('resolvedOutput: ', task.resolvedOutput)
//})

describe('Parallel', function() {
  it('should run two tasks in parallel', function(done) {
    const task1 = delayedWrite('foo', 1000)
    const task2 = delayedWrite('bar', 500)

    parallel(task1, task2)()
      .on('task.finish', function(task) {
//        console.log(task)
//        const data = this.output()
//
//        assert(data[0].toString() === 'bar')
//        assert(data[1].toString() === 'foo')
//
        done()
      })
  })

  it.skip('should parallel parallel', function(done) {
    const taskA1 = delayedResolve('A1', 100)
    const taskA2 = delayedResolve('A2', 500)
    const taskB1 = delayedResolve('B1', 200)
    const taskB2 = delayedResolve('B2', 400)

    const p1 = parallel(taskA1, taskA2)
    const p2 = parallel(taskB1, taskB2)

    parallel(p1, p2)()
      .on('close', function() {
        const data = this.output()

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
