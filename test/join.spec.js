'use strict'

const fs = require('fs')
const path = require('path')
const { assert } = require('chai')
const through = require('through2')
const intoStream = require('into-stream')
const split = require('split')

const Task = require('../lib/Task.js')
const join = require('../lib/Join.js')

const writeNumbers = (file) => Task({
  input: { value: ['1\n', '2\n'] },
  output: { file },
  name: 'Write numbers to file',
  skippable: false
}, ({ input }) =>
  intoStream(input)
    .pipe(fs.createWriteStream(file))
)

const sumNumbers = (numbers, sum) => Task({
  input: { file: numbers },
  output: { file: sum },
  name: 'Sum numbers in file',
  skippable: false
}, ({ input }) =>
  fs.createReadStream(input)
    .pipe(split())
    .pipe(through(function(chunk, enc, done) {
      if (!this._sum) {
        this._sum = 0
      }

      const line = chunk.toString()

      if (line !== '') {
        this._sum += parseInt(line)
      } else {
        this.push('' + this._sum)
      }

      done()
    }))
    .pipe(fs.createWriteStream(sum))
)

describe('Join', function() {
  it.skip('should join two tasks with stream, stream', function(done) {
    join(writeNumbers('numbers.1.txt'), sumNumbers('numbers.1.txt', 'sum.1.txt'))()
      .on('close', done)
      .on('error', done)
  })

  it.skip('should join two tasks with promise, promise', function(done) {

  })

  it.skip('should join two tasks with promise, stream', function(done) {

  })

  it.skip('should join two tasks with stream, promise', function(done) {

  })

  it.skip('should join joins', function(done) {
    const task1 = join(writeNumbers('numbers.2.txt'))
    const task2 = join(sumNumbers('numbers.2.txt', 'sum.2.txt'))

    join(task1, task2)()
      .on('close', done)
      .on('error', done)
  })
})
