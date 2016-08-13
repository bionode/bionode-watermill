'use strict'

const { assert } = require('chai')

const lifecycle = require('./lifecycle.js')

describe('lifecycle', function() {
  describe('required methods', function() {
    it('should have a create method', () => assert.hasOwnProperty('create'))
    it('should have a resumable method', () => assert.hasOwnProperty('resumable'))
    it('should have a resolveInput method', () => assert.hasOwnProperty('resolveInput'))
    it('should have a createOperation method', () => assert.hasOwnProperty('createOperation'))
    it('should have a settleOperation method', () => assert.hasOwnProperty('settleOperation'))
    it('should have a resolveOutput method', () => assert.hasOwnProperty('resolveOutput'))
    it('should have a validateOutput method', () => assert.hasOwnProperty('validateOutput'))
    it('should have a addDAGNode method', () => assert.hasOwnProperty('addDAGNode'))
  })

  it('each method should take a taskState and return a referentially dis-equal object', function(done) {
    for (const method in lifecycle) {
      // if (method === 'createOperation') break
      const taskState = {}
      const newTaskState = lifecycle[method](taskState)


      assert.notEqual(newTaskState, taskState, `lifecycle method ${method} did not return a new object`)
    }
    done()
  })
})
