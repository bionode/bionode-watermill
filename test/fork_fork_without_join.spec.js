'use strict'

const { task, fork, join, junction } = require('../')

const { assert } = require('chai')
const fs = require('fs-extra')

// tasks

const lsTask = (flags = '') => task({
  output: `*.txt`,
  params: { flags },
  name: `ls ${flags} 1`
}, ({ params }) => `ls ${params.flags} > files.txt`)

const ls = lsTask()
const lsAL = lsTask('-al')
const lsLAH = lsTask('-lah')

const lineCount = task({
  input: '*.txt',
  output: '*.count',
  name: 'Count lines from *.txt 1'
}, ({ input }) => `cat ${input} | wc -l > lines.count`)

const task0 = task({name: 'task0_1'}, () => `echo "something0"`)


describe('fork_fork_without_join', () => {
  it('should have a fork inside a fork without any wrapping', function (done) {
    const pipeline2 = join(
      task0,
      fork(
        fork(lsLAH, lsAL),
        ls
      ),
      lineCount
    )
    pipeline2().then((results) => {
      console.log('RESULTS: ', results)
      const obj = JSON.parse(fs.readFileSync('graphson.json', 'utf8'))
      // all these tests check if final task is the same in task definition
      // and trajectory
      assert.equal(results[0][0].tasks[1], results[0][0].context.trajectory[3])
      assert.equal(results[0][1].tasks[1], results[0][1].context.trajectory[3])
      assert.equal(results[1].tasks[1], results[1].context.trajectory[3])

      // This counts every vertex and edge created in this test and not only
      // this last pipeline2 graph
      assert.equal(obj.graph.vertices.length, 29) // 7 in this pipeline
      assert.equal(obj.graph.edges.length, 25) //6 in this pipeline

      done()
    })
  }).timeout(5000)
})
