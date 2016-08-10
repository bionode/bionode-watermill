'use strict'

const fs = require('fs')
const path = require('path')
const { assert } = require('chai')

let task, parallel, join

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
    join = waterwheelInstance.join
  })

  it('should run two tasks in parallel', function(done) {
    const task1 = delayedWrite('foo', 1000)
    const task2 = delayedWrite('bar', 500)

    parallel(task1, task2)().async.then((results) => {
      console.log('results: ', results)
      done()
    })
  })

  it('should parallel on join', function (done) {
    const taskA = delayedWrite('A', 100)
    const taskB = delayedWrite('B', 200)
    const taskC = delayedWrite('C', 150)
    const taskD = delayedWrite('D', 125)

    const pipeline = parallel(
      join(taskA, taskB),
      join(taskC, taskD)
    )

    pipeline().async.then((results) => {
      console.log('results: ', results)
      done()
    })
  })

  it('should parallel on parallel', function (done) {
    const taskA = delayedWrite('A', 100)
    const taskB = delayedWrite('B', 200)
    const taskC = delayedWrite('C', 150)
    const taskD = delayedWrite('D', 125)

    const pipeline = parallel(
      parallel(taskA, taskB),
      parallel(taskC, taskD)
    )

    pipeline().async.then((results) => {
      console.log('results: ', results)
      done()
    })
  })

  it('should parallel on join+parallel', function (done) {
    const taskA = delayedWrite('A', 100)
    const taskB = delayedWrite('B', 200)
    const taskC = delayedWrite('C', 150)
    const taskD = delayedWrite('D', 125)

    const pipeline = parallel(
      parallel(taskA, taskB),
      join(taskC, taskD)
    )

    pipeline().async.then((results) => {
      console.log('results: ', results)
      done()
    })
  })

  it('should arrive in proper order')
})
