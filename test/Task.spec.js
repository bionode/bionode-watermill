const { assert, expect, should } = require('chai')
const streamAssert = require('stream-assert')
const isStream = require('isstream')
const { isReadable, isWritable, isDuplex } = isStream

const Task = require('../lib/Task.js')

const noop = () => {}
const fs = require('fs')

describe('Task', function() {
  describe('checking parameters', function(done) {
    it.skip('should error if not provided with a function', function() {

    })
  })

  describe('when provided with a Promise', function() {
    describe('that resolves to a value', function() {
      it('should return a readable stream', function(done) {
        const task = Task(null, () => new Promise((resolve, reject) => resolve('datums')))

        assert.isOk(isReadable(task))
        assert.isNotOk(isWritable(task))

        task
          .on('data', noop)
          .on('end', () => done())
      })
    })
  })

  describe('when provided with a curried callback(err, data)', function() {
    it('should return a readable stream', function(done) {
      const task = Task(null, (props) => (cb) => cb(null, 'datums'))

      assert.isOk(isReadable(task))
      assert.isNotOk(isWritable(task))

      task
        .on('data', noop)
        .on('end', () => done())
    })
  })
})
