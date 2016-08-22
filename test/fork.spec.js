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

const lineCount = task({
  input: '*.txt',
  output: '*.count',
  name: 'Count lines from *.txt'
}, ({ input }) => `cat ${input} | wc -l > lines.count`)

describe('fork', function() {
  // Running this adds nodes to DAG breaking other tests
  // TODO setup/teadown
  it.skip('should return an array of results', function (done) {
    fork(ls, lsAL)().then((results) => {
      console.log('results: ', results)
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
      console.log('RESULTS: ', results)
      assert.equal(results[0].tasks[0], ls.info.uid)
      assert.equal(results[0].tasks[1], hash(lineCount.info.uid + 0))
      assert.equal(results[1].tasks[0], lsAL.info.uid)
      assert.equal(results[1].tasks[1], hash(lineCount.info.uid + 1))

      done()
    })
  })
})
