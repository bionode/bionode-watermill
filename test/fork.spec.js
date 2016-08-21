'use strict'

const { task, fork } = require('../')

const { assert } = require('chai')
const _ = require('lodash')
const path = require('path')
const twd = path.resolve(__dirname, 'files', 'fork')

const lsTask = (flags = '') => task({
  output: `${twd}/*.txt`,
  params: { flags }
}, ({ params }) => `ls ${params.flags} > ${twd}/${Date.now()}.txt`)

const ls = lsTask()
const lsAL = lsTask('-al')

describe('fork', function() {
  it ('should return an array of results', function (done) {
    fork(ls, lsAL)().then((results) => {
      console.log('results: ', results)
      assert.isOk(_.isArray(results))
      done()
    })
  })
})
