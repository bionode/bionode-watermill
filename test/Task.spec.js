'use strict'

const fs = require('fs')
const path = require('path')

const { assert, expect, should } = require('chai')
const streamAssert = require('stream-assert')
const isStream = require('isstream')
const { isReadable, isWritable, isDuplex } = isStream
const intoStream = require('into-stream')
const through = require('through2')
const duplexify = require('duplexify')
const split = require('split')

const { task, join, store } = require('../')()
const switchExt = require('../lib/utils/switch-ext.js')

const noop = () => {}

const twd = path.resolve(__dirname, 'files', 'task')

describe('Task', function() {
  describe('checking parameters', function(done) {
    it.skip('should error if not provided with a function', function() {

    })
  })

  describe('should interact with collection', function() {
    it('with no parameters', function(done) {
      const capitalize = task({
        input: '*.lowercase',
        output: '*.uppercase',
        name: `Capitalize *.lowercase -> *.uppercase`,
        dir: twd,
        resume: 'off'
      }, ({ input, dir }) =>
        fs.createReadStream(input)
          .pipe(through((chunk, enc, next) => next(null, chunk.toString().toUpperCase())))
          .pipe(fs.createWriteStream(switchExt(input, 'uppercase').replace(/source/, 'sink')))
      )

      join(capitalize)().on('task.finish', (results) => {
        const collection = store.getState().collection

        for (const uid of results.joined) {
          const task = store.getState().tasks[uid]

          assert.equal(collection.vertexValue(uid), task.resolvedOutput)
        }

        done()
      })
    })

    it('with parameters', function(done) {
      const capitalize = task({
        input: '*.lowercase',
        output: '*.uppercase',
        params: { foo: 'bar' },
        name: `Capitalize *.lowercase -> *.uppercase`,
        dir: twd,
        resume: 'off'
      }, ({ input, dir }) =>
        fs.createReadStream(input)
          .pipe(through((chunk, enc, next) => next(null, chunk.toString().toUpperCase())))
          .pipe(fs.createWriteStream(switchExt(input, 'uppercase').replace(/source/, 'sink')))
      )

      join(capitalize)().on('task.finish', (results) => {
        const collection = store.getState().collection

        for (const uid of results.joined) {
          const task = store.getState().tasks[uid]

          assert(collection.hasVertex(JSON.stringify(task.params)))
          assert(collection.hasEdge(JSON.stringify(task.params), uid))

          assert.equal(collection.vertexValue(uid), task.resolvedOutput)
        }

        done()
      })
    })
  })

  describe('when provided with a Promise', function() {
    describe('that resolves to a value', function() {
      it.skip('should return a readable stream', function(done) {
        const task = Task(null, () => new Promise((resolve, reject) => resolve('datums')))()

        assert.isOk(isReadable(task))
        // assert.isNotOk(isWritable(task))

        task
          .on('data', noop)
          .on('end', done)
      })
    })

    it.skip('should return a readable stream given a string', function(done) {
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


    it.skip('should return a readable stream given an object', function(done) {
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
    it.skip('should return a readable stream', function(done) {
      const task = Task(null, (props) => (cb) => cb(null, 'datums'))()

      assert.isOk(isReadable(task))
      // assert.isNotOk(isWritable(task))

      task
        .on('data', noop)
        .on('end', () => done())
    })
  })

  describe('when provided with a readable stream', function() {
    it.skip('should return a readable stream', function(done) {
      const task = Task(null, () => intoStream('unicorn'))()

      assert.isOk(isReadable(task))
      // assert.isNotOk(isWritable(task))

      task
        .on('data', noop)
        .on('end', () => done())
    })
  })

  describe('when provided with a writable stream', function() {
    it.skip('should return a writable stream', function(done) {
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
    it.skip('should return a duplex stream', function(done) {
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
