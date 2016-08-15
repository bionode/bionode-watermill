'use strict'

const { takeEvery } = require('redux-saga')
const { call, put, take, select } = require('redux-saga/effects')

function* lifecycle (action) {
  console.log('action: ', action)
  console.log('state from sagas: ', yield select())
}

function* rootSaga() {
  yield* takeEvery('tasks/create', lifecycle)
}

module.exports = rootSaga
