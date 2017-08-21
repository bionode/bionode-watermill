'use strict'

const { task, fork, join, junction } = require('../')

const { assert } = require('chai')
const fs = require('fs');

describe('fork_junction', () => {
  it('junction should work after a fork', (done) => {

    const task0 = task({name: 'task0_3"'}, () => `echo "something0_3"`)

    const task1 = task({name: 'task1_3'}, () => `echo "something1_3"`)

    const task2 = task({name: 'task2_3'}, () => `echo "something2_3"`)

    const task3 = task({name: 'task3_3'}, () => `echo "something3_3"`)

    const task4 = task({name: 'task4_3'}, () => `echo "something4_3"`)

    const task5 = task({name: 'task5_3'}, () => `echo "something5_3"`)

    const task6 = task({name: 'task6_3'}, () => `echo "something6_3"`)

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
      const obj = JSON.parse(fs.readFileSync('./graphson.json', 'utf8'));

      // should assure that junction node exists
      assert.isOk(results[0].context.trajectory[0])
      // should assure that both ends of the pipeline exist
      assert.isOk(results[0].context.trajectory[2])
      assert.isOk(results[1].context.trajectory[2])
      //checks if pipeline rendered the right number of vertices and edges
      assert.equal(obj.graph.vertices.length, 9)
      assert.equal(obj.graph.edges.length, 9)
      done()  // without done it is not really testing anything
    })
  }).timeout(5000)
})
