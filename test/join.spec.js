'use strict'

const { assert } = require('chai')

const { task, join } = require('../')

describe('Join', function() {
  it('Should pass context between two tasks', function(done) {
    const tasks = ['foo', 'bar'].map((val) =>
      task({
        params: { val },
        name: `Resolve ${val}`
      }, ({ params }) => new Promise(resolve => resolve(params.val))
      )
    )

    join(tasks[0], tasks[1])().then((results) => {
      assert.equal(results.context.trajectory[0], JSON.stringify({val: 'foo'}))
      assert.equal(results.context.trajectory[2], JSON.stringify({val: 'bar'}))
      assert.equal(results.type, 'join')

      done()
    })
  })
})
