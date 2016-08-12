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

/**
 * The most minimal `task` implementation and its responsibilties.
 *
 * An implementation of `bionode-waterwheel` shall expose a user facing function
 * called `task` which takes two parameters:
 * - the first, `props`, is a plain object describing properties of the task
 * - the second, `operationCreator`, is a function which receives one parameter:
 *   a plain object called `opsProps` which describes the properties of the
 *   operation. `opsProps` **must not be referentially equal to `props`**.
 * - returns a function `invocableTask` which will accept a callback `cb` and a
 *   context `ctx`. the callback shall be called on completion or error of the
 *   `operation`.
 */

const _ = require('lodash')
const noop = () => {}

/**
 * The task factory function.
 *
 * @param props {Object} properties of this task
 * @param operationCreator {Function} function to be called with props
 */
const task = (props, operationCreator) => {
  // Type checking on params
  if (!_.isPlainObject(props) || !_.isFunction(operationCreator)) {
    throw new Error('incorrect parameter types: task(props, operationCreator) where props is a plain object and operationCreator is a function')
  }

  // Returns an invocable task that can be passed a cb and optionally a ctx
  // By defauult cb is a noop and ctx is a object literal
  return (cb = () => {}, ctx = {}) => {
    const opsProps = Object.assign({}, props)
    const operation = operationCreator(opsProps)

    asyncDone(() => operation, (err, results) => {
      if (err) return cb(err)

      return cb(null, results)
    })
  }
}

// === TESTS ===
const { assert } = require('chai')
const asyncDone = require('async-done')

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
      assert(_.isFunction(task({}, noop)), 'task did not return a function')
    })

    it('calls operationCreator with an opsProps object', function(done) {
      task({}, (opsProps) => {
        assert.isOk(_.isPlainObject(opsProps), 'opsProps was not a plain object')
        done()
      })()
    })

    it('opsProps should not referentially equal props', function(done) {
      const props = {}
      task(props, (opsProps) => {
        assert.isOk(props !== opsProps)
        done()
      })()
    })

    it('should eventually call callback', function(done) {
      task(
        { payload: 'foo' },
        ({ payload }) => new Promise(resolve => resolve(payload.toUpperCase()))
      )( (err, data) => done() )
    })
  })
})


