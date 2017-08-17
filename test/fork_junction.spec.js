'use strict'

const { task, fork, join, junction } = require('../')

const { should } = require('chai')

const task0 = task({name: 'task0_2"'}, () => `echo "something0"`)

const task1 = task({name: 'task1'}, () => `echo "something1"`)

const task2 = task({name: 'task2'}, () => `echo "something2"`)

const task3 = task({name: 'task3'}, () => `echo "something3"`)

const task4 = task({name: 'task4'}, () => `echo "something4"`)

const task5 = task({name: 'task5'}, () => `echo "something5"`)

const task6 = task({name: 'task6'}, () => `echo "something6"`)

describe('fork_junction', function() {
  it('junction should work after a fork', function() {
    const pipeline = join(
      task0,
      fork(
        join(
          task1,
          junction(
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
      // should assure that junction node exists
      should.exist(results[0].context.trajectory[0])
      // should assure that both ends of the pipeline exist
      should.exist(results[0].context.trajectory[2])
      should.exist(results[1].context.trajectory[2])

      //done()
    })
  })
})
