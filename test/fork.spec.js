'use strict'

const { task, fork, join, junction } = require('../')

const { assert } = require('chai')
const _ = require('lodash')

const lsTask = (flags = '') => task({
  output: `*.txt`,
  params: { flags },
  name: `ls ${flags}`
}, ({ params }) => `ls ${params.flags} > files.txt`)

const ls = lsTask()
const lsAL = lsTask('-al')

const task1 = task({name: 'task1'}, () => `echo "something1"`)

const task2 = task({name: 'task2'}, () => `echo "something2"`)

describe('fork', () => {
  // Running this adds nodes to DAG breaking other tests
  // TODO setup/teadown
  it.skip('should return an array of results', (done) => {
    fork(ls, lsAL)().then((results) => {
      //console.log('results: ', results)
      assert.isOk(_.isArray(results))
      done()
    })
  })

  it ('should have info.type be "fork"', (done) => {
    const forked = fork(ls, lsAL)

    assert.equal(forked.info.type, 'fork')
    done()
  })

  it ('should work with join', (done) => {
    const forked = fork(ls, lsAL)
    const fork_join = join(task1, forked, task2)
    fork_join().then((results) => {
      console.log('RESULTS: ', results)
      assert.equal(results[0].tasks[0], ls.info.uid)
      // newUid of task is inaccessible because context doesn't exist prior
      // to task execution
      // instead the resulting new task id (newUid from task.js) is added to
      // task.context.trajectory
      assert.equal(results[0].tasks[1], results[0].context.trajectory[3])
      assert.equal(results[1].tasks[0], lsAL.info.uid)
      // newUid of task is unnaccessible because context doesn't exist prior'
      //' to task execution
      assert.equal(results[1].tasks[1], results[1].context.trajectory[3])
    })
    done()
  })
})
