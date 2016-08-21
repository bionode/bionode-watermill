'use strict'

const { takeEvery } = require('redux-saga')

const { CREATE_TASK } = require('../reducers/tasks.js')
const lifecycle = require('./lifecycle.js')

function* rootSaga() {
  yield* takeEvery(CREATE_TASK, lifecycle)
}

module.exports = rootSaga

