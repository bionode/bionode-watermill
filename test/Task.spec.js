'use strict'

/**
 * The task function takes in `props` and an `operationCreator`. A series of
 * transformations will be applied to the props through the lifecycle methods
 * which are responsible for receiving props and mutating them. The quick way is
 * to modify the same `props` object over a closure; the safe way is to use a pure
 * function that takes an immutable `props` and returns an `opsProps` (like a
 * redux reducer). However, this detail is internal to the abstract `task`
 * implementation for our purposes (the API definition and test suite). `task`
 * returns a function which takes a `ctx`. The `ctx` is used to store the
 * position of this task in the DAG (`ctx.trajectory`, where `trajectory` is an
 * array of vertex ids).
 *
 * Ultimately, an invocation to `task` shall return a function which returns
 * a function, stream, promise, observable, callback - anything that can be
 * handled by async-done. In this way - think of the item that `task()()`
 * produces to be a valid function for a Gulp task.
 */

const { assert } = require('chai')
const _ = require('lodash')

const { task } = require('../')

describe('task', function() {
  describe('parameter type checks', function() {
    it('should throw if props is not a plain object', function() {
      assert.throws(() => task('foo', noop), Error)
    })

    it('should throw if operationCreator is not a function', function() {
      assert.throws(() => task({}, 'foo'), Error)
    })
  })

  describe('invocable task', function() {
    it('should be a function', function() {
      assert(_.isFunction(task({}, _.noop)), 'task did not return a function')
    })

    it('calls operationCreator with an opsProps object', function(done) {
      task({ resume: 'off' }, (opsProps) => {
        assert.isOk(_.isPlainObject(opsProps), 'opsProps was not a plain object')
        done()
      })()
    })

    it('opsProps should not referentially equal props', function(done) {
      const props = { resume: 'off' }
      task(props, (opsProps) => {
        assert.isOk(props !== opsProps)
        done()
      })()
    })

    it('should eventually call callback with resume on', function(done) {
      task(
        { payload: 'foo' },
        ({ payload }) => new Promise(resolve => resolve(payload.toUpperCase()))
      )( (err, data) => done() )
    })

    it('should eventually call callback with resume off', function(done) {
      task(
        { payload: 'foo', resume: 'off' },
        ({ payload }) => new Promise(resolve => resolve(payload.toUpperCase()))
      )( (err, data) => done() )
    })
  })
})


