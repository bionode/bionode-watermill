'use strict'

const { task, fork, join, junction } = require('../')

const { assert } = require('chai')
const _ = require('lodash')
const path = require('path')
const twd = path.resolve(__dirname, 'files', 'fork')
const hash = require('../lib/utils/hash.js')

const lsTask = (flags = '') => task({
  output: `*.txt`,
  params: { flags },
  name: `ls ${flags}`
}, ({ params }) => `ls ${params.flags} > files.txt`)

const ls = lsTask()
const lsAL = lsTask('-al')
const lsLAH = lsTask('-lah')

const lineCount = task({
  input: '*.txt',
  output: '*.count',
  name: 'Count lines from *.txt'
}, ({ input }) => `cat ${input} | wc -l > lines.count`)

const task0 = task({name: 'task0'}, () => `echo "something0"`)

describe('fork', function() {
  // Running this adds nodes to DAG breaking other tests
  // TODO setup/teadown
  it.skip('should return an array of results', function (done) {
    fork(ls, lsAL)().then((results) => {
      //console.log('results: ', results)
      assert.isOk(_.isArray(results))
      done()
    })
  })

  it ('should have info.type be "fork"', function() {
    const forked = fork(ls, lsAL)

    assert.equal(forked.info.type, 'fork')
  })

  it('should work with join', function(done) {
    const pipeline = join(fork(ls, lsAL), lineCount)

    pipeline().then((results) => {
      //console.log('RESULTS: ', results)
      assert.equal(results[0].tasks[0], ls.info.uid)
      // newUid of task is inaccessible because context doesn't exist prior
      // to task execution
      // instead the resulting new task id (newUid from task.js) is added to
      // task.context.trajectory
      assert.equal(results[0].tasks[1], results[0].context.trajectory[2])
      assert.equal(results[1].tasks[0], lsAL.info.uid)
      // newUid of task is unnaccessible because context doesn't exist prior to task execution
      assert.equal(results[1].tasks[1], results[1].context.trajectory[2])

      done()
    })
  })

  it('should have a fork inside a fork without any wrapping', function(done) {
    const pipeline2 = join(
      task0,
      fork(
        fork(lsLAH, lsAL),
        ls
      ),
      lineCount
    )
    pipeline2().then((results) => {
      //console.log('RESULTS: ',results)
      // all these tests check if final task is the same in task definition
      // and trajectory
      assert.equal(results[0][0].tasks[1], results[0][0].context.trajectory[3])
      assert.equal(results[0][1].tasks[1], results[0][1].context.trajectory[3])
      assert.equal(results[1].tasks[1], results[1].context.trajectory[3])

      done()
    })
  }).timeout(5000)
})
