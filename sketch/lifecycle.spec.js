'use strict'

const _ = require('lodash')
const { assert } = require('chai')

const lifecycle = require('./lifecycle')

describe('lifecycle', function() {
  describe('required methods', function() {
    it('should have a create method', () => assert.hasOwnProperty('create'))
    it('should have a resumable method', () => assert.hasOwnProperty('resumable'))
    it('should have a resolveInput method', () => assert.hasOwnProperty('resolveInput'))
    it('should have a createOperation method', () => assert.hasOwnProperty('createOperation'))
    it('should have a settleOperation method', () => assert.hasOwnProperty('settleOperation'))
    it('should have a resolveOutput method', () => assert.hasOwnProperty('resolveOutput'))
    it('should have a validateOutput method', () => assert.hasOwnProperty('validateOutput'))
    it('should have a postValidation method', () => assert.hasOwnProperty('postValidation'))
  })

  it('each method should take a taskState and return a referentially dis-equal object', function() {
    for (const method in lifecycle) {
      const taskState = {}
      const args = [ taskState ]
      if (method === 'createOperation') args.push(_.noop)
      const newTaskState = lifecycle[method].apply(null, args)

      assert.notEqual(newTaskState, taskState, `lifecycle method ${method} did not return a new object`)
    }
  })
})
