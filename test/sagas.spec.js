'use strict'

const { assert } = require('chai')
const fs = require('fs')
const through = require('through2')

const { select, call } = require('redux-saga/effects')

const { defaultTask } = require('../lib/constants/default-task-state.js')
const lifecycle = require('../lib/sagas/lifecycle.js')
const {
  generateUid,
  checkResumable
} = require('../lib/lifecycle')
const { createTask, CREATE_TASK } = require('../lib/reducers/tasks.js')
const { selectTask } = require('../lib/sagas/selectors.js')

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

  const generator = lifecycle(createTask({
    uid,
    hashes,
    props,
    operationCreator,
    taskResolve: resolve,
    taskReject: reject
  }))

  let next = generator.next()
  it.skip('Should select for the task', function() {
    assert.isOk(next.value.SELECT.function === select(selectTask(uid)).SELECT.function)
  })

  const taskState = Object.assign({}, defaultTask, props, { uid, hashes })
  it.skip('Should first check if resumable', function() {
    next = generator.next(taskState)
    assert.deepEqual(next.value, call(checkResumable, taskState))
  })
})
