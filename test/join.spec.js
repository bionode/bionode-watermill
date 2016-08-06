'use strict'

const fs = require('fs')
const path = require('path')
const { assert } = require('chai')
const through = require('through2')
const intoStream = require('into-stream')
const split = require('split')

const { task, join } = require('../')

const writeNumbers = (file) => task({
  input: null,
  output: file,
  name: `Write numbers to ${file}`,
  params: { numbers: ['1\n', '2\n'] },
  resume: 'off'
}, ({ params }) => intoStream(params.numbers).pipe(fs.createWriteStream(file)) )

const sumNumbers = (numbers, sum) => task({
  input: numbers,
  output: sum,
  name: `Sum numbers in ${numbers}`,
  resume: 'off'
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
  // Skipping b/c breaks when more than one are ran
  it.skip('should join two tasks with stream, stream', function(done) {
    join(writeNumbers('numbers.1.txt'), sumNumbers('numbers.1.txt', 'sum.1.txt'))()
      .on('task.finish', (task) => {
        done()
      })
  })

  it.skip('should join two tasks with promise, promise', function(done) {

  })

  it.skip('should join two tasks with promise, stream', function(done) {

  })

  it.skip('should join two tasks with stream, promise', function(done) {

  })

  it('should join joins', function(done) {
    const task1 = join(writeNumbers('numbers.2.txt'))
    // const task2 = join(sumNumbers('numbers.2.txt', 'sum.2.txt'))
    const task2 = sumNumbers('numbers.2.txt', 'sum.2.txt')

    join(task1, task2)()
      .on('task.finish', () => done())
  })
})
