const { assert, expect, should } = require('chai')
const streamAssert = require('stream-assert')
const isStream = require('isstream')
const { isReadable, isWritable, isDuplex } = isStream
const intoStream = require('into-stream')
const through = require('through2')
const duplexify = require('duplexify')
const split = require('split')

const Task = require('../lib/Task.js')
const join = require('../lib/Join.js')

const noop = () => {}
const fs = require('fs')
const path = require('path')

describe('Task', function() {
  describe('checking parameters', function(done) {
    it.skip('should error if not provided with a function', function() {

    })
  })

  describe('when provided with a Promise', function() {
    describe('that resolves to a value', function() {
      it('should return a readable stream', function(done) {
        const task = Task(null, () => new Promise((resolve, reject) => resolve('datums')))()

        assert.isOk(isReadable(task))
        assert.isNotOk(isWritable(task))

        task
          .on('data', noop)
          .on('end', done)
      })
    })

    it('should return a readable stream given a string', function(done) {
      const task = Task({
        input: {
          one: { value: 5 },
          two: { value: 2 }
        },
        name: 'testing promise stream with string'
      }, () => new Promise((resolve, reject) => resolve('42')))
      

      task()
        .on('data', (data) => {
          assert(data.toString() === '42')
          
          done() 
        })
        // .on('end', done)
        .on('error', done)
    })


    it('should return a readable stream given an object', function(done) {
      const task = Task({
        input: {
          one: { value: 5 },
          two: { value: 2 }
        },
        output: { stream: 'object' },
        name: 'testing promise stream with object' 
      }, () => new Promise((resolve, reject) => resolve({ foo: 'bar' })))
      

      task()
        .on('data', (data) => {
          assert.deepEqual(data, { foo: 'bar' })
          
          done()  
        })
        .on('error', done)
    }) 
  })

  describe('when provided with a curried callback(err, data)', function() {
    it('should return a readable stream', function(done) {
      const task = Task(null, (props) => (cb) => cb(null, 'datums'))()

      assert.isOk(isReadable(task))
      assert.isNotOk(isWritable(task))

      task
        .on('data', noop)
        .on('end', () => done())
    })
  })

  describe('when provided with a readable stream', function() {
    it('should return a readable stream', function(done) {
      const task = Task(null, () => intoStream('unicorn'))()

      assert.isOk(isReadable(task))
      assert.isNotOk(isWritable(task))

      task
        .on('data', noop)
        .on('end', () => done())
    })
  })

  describe('when provided with a writable stream', function() {
    it('should return a writable stream', function(done) {
      const task = Task(null, () => fs.createWriteStream(path.resolve(__dirname, 'writable.log')))()

      task.write('one\n')
      task.write('two\n')
      task.end()

      assert.isOk(isWritable(task))
      assert.isNotOk(isReadable(task))

      task
        .on('finish', () => done())
    })
  })

  describe('when provided with a duplex stream', function() {
    it('should return a duplex stream', function(done) {
      // Simple pass-through duplex stream
      const task = Task(null, () => through())()

      assert.isOk(isDuplex(task))

      task.write('datum')
      task.end()

      task
        .on('data', noop)
        //.on('end', () => done())
        .on('finish', () => done())
    })
  })
})

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
      .on('task.done', done)
      .on('error', done)  
  })

  it.skip('should join two tasks with promise, promise', function(done) {

  })

  it.skip('should join two tasks with promise, stream', function(done) {
    
  })

  it.skip('should join two tasks with stream, promise', function(done) {

  })
})
