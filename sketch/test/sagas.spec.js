'use strict'

const { assert } = require('chai')
const fs = require('fs')
const through = require('through2')

const { take, select } = require('redux-saga/effects')

const { lifecycle } = require('../lib/sagas')
const { generateUid } = require('../lib/lifecycle')
const { createTask, CREATE_TASK } = require('../lib/reducers/tasks.js')

describe('Task lifecycle saga', function() {
  const props = {
    input: '*.source',
    output: '*.sink'
  }
  const operationCreator = ({ input }) =>
    fs.createReadStream(input)
      .pipe(through((chunk, enc, cb) => cb(null, chunk.toString().toUpperCase())))
      .pipe(fs.createWriteStream(input.replace(/\.source$/, 'sink')))

  // This is the portion handled when you first call task()
  const { uid, hashes } = generateUid(props)
  const resolve = () => console.log('Placeholder promise resolve')
  const reject = () => console.log('Placeholder promise reject')

  const iterator = lifecycle(createTask({
    uid,
    hashes,
    props,
    operationCreator,
    taskResolve: resolve,
    taskReject: reject
  }))

  it('then it should select', function() {
    assert.deepEqual(
      iterator.next().value,
      select()
    )
  })
})
