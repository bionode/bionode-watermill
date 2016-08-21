'use strict'

const _ = require('lodash')
const { assert } = require('chai')

const statusTypes = require('../lib/constants/task-status-types.js')
const defaultConfig = require('../lib/constants/default-config-state.js')
const lifecycle = require('../lib/lifecycle')

describe('lifecycle', function() {
  describe('required methods', function() {
    [
      // 'create',
      // 'resumable',
      'resolveInput',
      'resolveOutput',
      'createOperation',
      'settleOperation',
      'resolveOutput',
      'validateOutput',
      'postValidation'
    ].map((method) => {
      it(`should have a ${method} method`, function() {
        assert.isOk(lifecycle.hasOwnProperty(method))
      })
    })
  })

  it.skip('each method should take a taskState and return a referentially dis-equal object', function() {
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
        it.skip('should have proper types', function(done) {
          lifecycle.create({}).then((results) => {
            assertProperTypes(mappings, results)
            done()
          })
        })
      })

      describe('when props has threads', function() {
        it.skip('should use passed in value', function() {
          const threads = 2
          lifecycle.create({ threads }).then((results) => {
            assert.equal(initialTaskState.threads, threads)

            assertProperTypes(mappings, initialTaskState)

            done()
          })
        })

        it.skip('should throw if not a number', function(done) {
          const threads = '2'

          lifecycle.create({ threads }).catch((err) => {
            // TODO assert err is 'Bad type for threads, got string but expected number'
            done()
          })
        })
      })
    })

    describe('resolveInput', function() {
      it.skip('should return an object with key resolvedInput', function(done) {
        lifecycle.create().then(taskState => {
          lifecycle.resolveInput(taskState).then((results) => {
            assert.isOk(results.hasOwnProperty('resolvedInput'))
            done()
          })
        })
      })
    })

    describe('createOperation', function() {
      it.skip('should return an object with keys operation and operationString', function(done) {
        lifecycle.create().then(taskState => {
          lifecycle.createOperation(taskState, () => {}).then((results) => {
            assert.isOk(results.hasOwnProperty('operation'))
            assert.isOk(results.hasOwnProperty('operationString'))
            done()
          })
        })
      })
    })

    describe('settleOperation', function() {
      it.skip('should call callback', function(done) {
        lifecycle.settleOperation(new Promise(resolve => resolve('foo'))).then((results) => {
          done()
        })
      })
    })

    describe('resolveOutput', function() {
      it.skip('should return an object with key resolvedOutput', function(done) {
          lifecycle.create().then(taskState => {
          lifecycle.resolveOutput(taskState).then((results) => {
            assert.isOk(results.hasOwnProperty('resolvedOutput'))
          })
          done()
        })
      })
    })

    describe('validateOutput', function() {
      it.skip('should return an object with keys validations and validated', function(done) {
        lifecycle.create().then(taskState => {
          lifecycle.validateOutput(taskState).then((results) => {
            assert.isOk(result.hasOwnProperty('validations'))
            assert.isOk(_.isPlainObject(result.validations))
            assert.isOk(result.hasOwnProperty('validated'))
            assert.isBoolean(result.validated)
          })

          done()
        })
      })
    })

    describe('postValidation', function() {
      it.skip('should return an object with key postValidations', function(done) {
        lifecycle.create().then(taskState => {
          lifecycle.postValidation(taskState).then((results) => {
            assert.isOk(result.hasOwnProperty('postValidation'))
            assert.isOk(_.isPlainObject(result.postValidation))
          })

          done()
        })
      })
    })
  })
})
