'use strict'

const { task, fork, join, junction } = require('../')

const { assert } = require('chai')

const task0 = task({name: 'task0_3"'}, () => `echo "something0"`)

const task1 = task({name: 'task1_3'}, () => `echo "something1"`)

const task2 = task({name: 'task2_3'}, () => `echo "something2"`)

const task3 = task({name: 'task3_3'}, () => `echo "something3"`)

const task4 = task({name: 'task4_3'}, () => `echo "something4"`)

const task5 = task({name: 'task5_3'}, () => `echo "something5"`)

const task6 = task({name: 'task6_3'}, () => `echo "something6"`)

describe('fork_fork', function() {
  it('fork should work after a fork, wrapped in a join', function() {
    const pipeline = join(
      task0,
      fork(
        join(
          task1,
          fork(
            task4,
            task5
          ),
          task6
        ),
        task2
      ),
      task3
    )
    pipeline().then((results) => {
      console.log('RESULTS: ',results)

      // should assure that both ends of the pipeline exist
      assert.equal(results[1].tasks[1],results[1].context.trajectory[2])
    })
  })
})
