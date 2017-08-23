'use strict'

const { task, fork, join, junction } = require('../')

const { assert } = require('chai')
const fs = require('fs-extra')

// Define tasks

const lsTask = (flags = '') => task({
  output: `*.txt`,
  params: { flags },
  name: `ls ${flags} 2`
}, ({ params }) => `ls ${params.flags} > files.txt`)

const ls = lsTask()
const lsAL = lsTask('-al')
const lsLAH = lsTask('-lah')

const lineCount = task({
  input: '*.txt',
  output: '*.count',
  name: 'Count lines from *.txt 2'
}, ({ input }) => `cat ${input} | wc -l > lines.count`)

const task0 = task({name: 'task0_2'}, () => `echo "something0"`)

const task1 = task({name: 'task1_2'}, () => `echo "something1"`)

const task2 = task({name: 'task2_2'}, () => `echo "something2"`)


describe('fork_junction', () => {
  it('junction should work after a fork', (done) => {
    const pipeline3 = join(
      task0,
      fork(
        join(
          task1,
          junction(
            ls,
            lsAL
          ),
          task2
        ),
        lsLAH
      ),
      lineCount
    )
    pipeline3().then((results) => {
      console.log('RESULTS: ',results)
      const obj = JSON.parse(fs.readFileSync('graphson.json', 'utf8'))

      // should assure that junction node exists
      assert.isOk(results[0].context.trajectory[0])
      // should assure that both ends of the pipeline exist
      assert.isOk(results[0].context.trajectory[2])
      assert.isOk(results[1].context.trajectory[2])
      //checks if pipeline rendered the right number of vertices and edges
      assert.equal(obj.graph.vertices.length, 39)
      assert.equal(obj.graph.edges.length, 35)
      done()  // without done it is not really testing anything
    })
  }).timeout(5000)
})
