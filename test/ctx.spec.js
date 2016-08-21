'use strict'

const { assert } = require('chai')

// const { defaultContext } = require('../lib/constants/default-task-state.js')
const { mergeCtx } = require('../lib/ctx')

describe('Context', function() {
  describe('mergeCtx', function() {
    it('should merge two ctx as expected', function() {
      const currentCtx = {
        uid: 'one',
        trajectory: []
      }
      const results = {
        uid: 'two',
        params: { val: 'foo' },
        context: {
          trajectory: []
        }
      }

      const mergedCtx = mergeCtx('type')(currentCtx, results)

      assert.deepEqual(mergedCtx, {
        trajectory: [
          JSON.stringify(results.params),
          'two'
        ],
        type: 'type'
      })
    })
  })
})
