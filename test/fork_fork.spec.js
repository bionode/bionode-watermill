'use strict'

const { task, fork, join, junction } = require('../')

const { assert } = require('chai')
const fs = require('fs-extra')
console.log("teste fork fork")

// Define tasks

const lsTask = (flags = '') => task({
  output: `*.txt`,
  params: { flags },
  name: `ls ${flags} 3`
}, ({ params }) => `ls ${params.flags} > files.txt`)

const ls = lsTask()
const lsAL = lsTask('-al')
const lsLAH = lsTask('-lah')

const lineCount = task({
  input: '*.txt',
  output: '*.count',
  name: 'Count lines from *.txt 3'
}, ({ input }) => `cat ${input} | wc -l > lines.count`)

const task0 = task({name: 'task0_3'}, () => `echo "something0"`)

// const task1 = task({name: 'task1_3'}, () => `echo "something1"`)
//
// const task2 = task({name: 'task2_3'}, () => `echo "something2"`)
//
// const task3 = task({name: 'task3_3'}, () => `echo "something3"`)

const task4 = task({name: 'task4_3'}, () => `echo "something4"`)

const task5 = task({name: 'task5_3'}, () => `echo "something5"`)

const task6 = task({name: 'task6_3'}, () => `echo "something6"`)

describe('fork_fork', () => {
  it('fork should work after a fork, wrapped in a join', (done) => {
    const pipeline4 = join(
      task0,
      fork(
        join(
          task4,
          fork(
            ls,
            lsAL
          ),
          task6
        ),
        lsLAH
      ),
      task5
    )

    pipeline4().then((results) => {
      console.log('RESULTS test: ',results)
      const obj = JSON.parse(fs.readFileSync('graphson.json', 'utf8'))

      // should assure that both ends of the pipeline exist
      //assert.equal(results[1].tasks[1], results[1].context.trajectory[2])
      //checks if pipeline rendered the right number of vertices and edges
      assert.equal(obj.graph.vertices.length, 17)
      assert.equal(obj.graph.edges.length, 14)
      done() // without done it is not really testing anything
    })
  }).timeout(5000)
})
