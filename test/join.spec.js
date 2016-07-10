'use strict'

const fs = require('fs')
const path = require('path')
const { assert } = require('chai')
const through = require('through2')
const intoStream = require('into-stream')
const split = require('split')

const Task = require('../lib/Task.js')
const join = require('../lib/Join.js')


describe('Join', function() {
  it('should join two tasks with stream, stream', function(done) {
    const task1 = Task({
      input: { value: ['1\n', '2\n'] },
      output: { file: 'foo.txt' },
      name: 'Write numbers to file'
    }, ({ input }) =>
      intoStream(input)
        .pipe(fs.createWriteStream('foo.txt'))
    )

    const task2 = Task({
      input: { file: 'foo.txt' },
      output: { file: 'sum.txt' },
      name: 'Sum numbers in file'
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
        .pipe(fs.createWriteStream('sum.txt'))
    )

    join(task1, task2)()
      .on('close', done)
      .on('error', done)  
  })

  it.skip('should join two tasks with promise, promise', function(done) {

  })

  it.skip('should join two tasks with promise, stream', function(done) {
    
  })

  it.skip('should join two tasks with stream, promise', function(done) {

  })
})
