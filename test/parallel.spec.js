'use strict'

const { assert } = require('chai')

const { task, parallel } = require('../')

describe('parallel', function() {
  it('should produce a concatenated trajectory', function(done) {
    const tasks = ['foo', 'bar'].map((val) =>
      task({
        params: { val }
      }, ({ params }) => new Promise(resolve => resolve(params.val))
      )
    )

    parallel.apply(null, tasks)().then((results) => {
      assert.equal(results.trajectory[0], JSON.stringify({val: 'foo'}))
      assert.equal(results.trajectory[2], JSON.stringify({val: 'bar'}))
      assert.equal(results.type, 'parallel')

      done()
    })
  })
})
