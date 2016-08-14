'use strict'

const { assert } = require('chai')

const task = require('./task.js')
const join = require('./orchestrators/join.js')

describe('Join', function() {
  it('Should pass context between two tasks', function(done) {
    const tasks = ['foo', 'bar'].map((val) =>
      task({
        params: { val }
      }, ({ params }) => new Promise(resolve => resolve(params.val))
      )
    )

    join(tasks[0], tasks[1])().then((results) => {
      assert.equal(results.trajectory[0], JSON.stringify({val: 'foo'}))
      assert.equal(results.trajectory[2], JSON.stringify({val: 'bar'}))

      done()
    })
  })
})
