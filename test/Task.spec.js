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
      it('should return a readable stream', function() {
        const task = Task(null, () => new Promise((resolve, reject) => resolve('datums')))

        assert.isOk(isReadable(task))
        assert.isNotOk(isWritable(task))
      })
    })
  })
})
