'use strict'

const { assert } = require('chai')

const { defaultCtx, mergeCtx } = require('../lib/ctx')

describe('Context', function() {
  describe('mergeCtx', function() {
    it('should merge two ctx as expected', function() {
      const ctx1 = {
        uid: 'one',
        trajectory: []
      }
      const ctx2 = {
        uid: 'two',
        params: { val: 'foo' },
        trajectory: []
      }

      const mergedCtx = mergeCtx('type')(ctx1, ctx2)

      assert.deepEqual(mergedCtx, {
        trajectory: [
          JSON.stringify(ctx2.params),
          'two'
        ],
        type: 'type'
      })
    })
  })
})
