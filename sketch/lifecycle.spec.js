'use strict'

const _ = require('lodash')
const { assert } = require('chai')

const statusTypes = require('./constants/task-status-types.js')
const defaultConfig = require('./constants/default-config-state.js')
const lifecycle = require('./lifecycle')

describe('lifecycle', function() {
  describe('required methods', function() {
    it('should have a create method', () => assert.isOk(lifecycle.hasOwnProperty('create')))
    it('should have a resumable method', () => assert.isOk(lifecycle.hasOwnProperty('resumable')))
    it('should have a resolveInput method', () => assert.isOk(lifecycle.hasOwnProperty('resolveInput')))
    it('should have a createOperation method', () => assert.isOk(lifecycle.hasOwnProperty('createOperation')))
    it('should have a settleOperation method', () => assert.isOk(lifecycle.hasOwnProperty('settleOperation')))
    it('should have a resolveOutput method', () => assert.isOk(lifecycle.hasOwnProperty('resolveOutput')))
    it('should have a validateOutput method', () => assert.isOk(lifecycle.hasOwnProperty('validateOutput')))
    it('should have a postValidation method', () => assert.isOk(lifecycle.hasOwnProperty('postValidation')))
  })

  it('each method should take a taskState and return a referentially dis-equal object', function() {
    for (const method in lifecycle) {
      const taskState = {}
      const args = [ taskState ]
      if (method === 'createOperation') args.push(_.noop)
      if (method === 'settleOperation') args.push(_.noop)
      const newTaskState = lifecycle[method].apply(null, args)

      assert.notEqual(newTaskState, taskState, `lifecycle method ${method} did not return a new object`)
    }
  })

  describe('methods', function() {
    describe('create', function() {

      const mappings = {
        threads: { it: 'should be a number', fn: assert.isNumber },
        container: { it: 'should be null', fn: assert.isNull  },
        resume: { it: 'should be on', fn: (val) => assert.equal('on', val) },
        uid: { it: 'should be a string', fn: assert.isString },
        hashes: { it: 'should be an object with strings input, output, params ', fn: (val) => {
          assert.isOk(val.hasOwnProperty('input'))
          assert.isString(val['input'])
          assert.isOk(val.hasOwnProperty('output'))
          assert.isString(val['output'])
          assert.isOk(val.hasOwnProperty('params'))
          assert.isString(val['params'])
        } },
        name: { it: 'should be "Unnamed Task"', fn: (val) => assert.equal('Unnamed Task', val) },
        dir: { it: 'should be null', fn: assert.isNull },
        input: { it: 'should be null', fn: assert.isNull },
        output: { it: 'should be null', fn: assert.isNull },
        operation: { it: 'should be null', fn: assert.isNull },
        status: { it: 'should be STATUS_CREATED', fn: (val) => assert.equal(statusTypes.STATUS_CREATED, val) },
        created: { it: 'should be a number', fn: assert.isNumber },
        validated: { it: 'should be false', fn: assert.isFalse },
        params: { it: 'should be {}', fn: (val) => assert.deepEqual({}, val) },
        trajectory: { it: 'should be []', fn: (val) => assert.deepEqual([], val) }
      }

      const assertProperTypes = (mappings, taskState) => {
        for (const key in mappings) {
          const val = mappings[key]

          val.fn.call(null, taskState[key])
        }
      }

      describe('when props is {}', function() {
        it('should have proper types', function() {
          const initialTaskState = lifecycle.create({})

          assertProperTypes(mappings, initialTaskState)
        })
      })

      describe('when props has threads', function() {
        it('should use passed in value', function() {
          const threads = 2
          const initialTaskState = lifecycle.create({ threads })

          assert.equal(initialTaskState.threads, threads)

          assertProperTypes(mappings, initialTaskState)
        })

        it('should throw if not a number', function() {
          const threads = '2'

          assert.throws(() => lifecycle.create({ threads }), 'Bad type for threads, got string but expected number')
        })
      })
    })

    describe('resumable', function() {
      it('should return a boolean', function() {
        const taskState = lifecycle.create({})
        const result = lifecycle.resumable(taskState)

        assert.isBoolean(result)
      })
    })

    describe('resolveInput', function() {
      it('should return an object with key resolvedInput', function() {
        const taskState = lifecycle.create()

        const result = lifecycle.resolveInput(taskState)
        assert.isOk(result.hasOwnProperty('resolvedInput'))
      })
    })

    describe('createOperation', function() {
      it('should return an object with keys operation and operationString', function() {
        const taskState = lifecycle.create()
        const result = lifecycle.createOperation(taskState, () => {})

        assert.isOk(result.hasOwnProperty('operation'))
        assert.isOk(result.hasOwnProperty('operationString'))
      })
    })

    describe('settleOperation', function() {
      it('should call callback', function(done) {
        lifecycle.settleOperation(
          () => new Promise(resolve => resolve('foo')),
          (err, data) => done()
        )
      })
    })

    describe('resolveOutput', function() {
      it('should return an object with key resolvedOutput', function() {
        const taskState = lifecycle.create()

        const result = lifecycle.resolveOutput(taskState)
        assert.isOk(result.hasOwnProperty('resolvedOutput'))
      })
    })

    describe('validateOutput', function() {
      it('should return an object with keys validations and validated', function() {
        const taskState = lifecycle.create()
        const result = lifecycle.validateOutput(taskState)

        assert.isOk(result.hasOwnProperty('validations'))
        assert.isOk(_.isPlainObject(result.validations))
        assert.isOk(result.hasOwnProperty('validated'))
        assert.isBoolean(result.validated)
      })
    })

    describe('postValidation', function() {
      it('should return an object with key postValidations', function() {
        const taskState = lifecycle.create()
        const result = lifecycle.postValidation(taskState)

        assert.isOk(result.hasOwnProperty('postValidation'))
        assert.isOk(_.isPlainObject(result.postValidation))
      })
    })
  })
})
