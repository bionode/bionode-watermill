'use strict'

const { addOutput } = require('../reducers/collection.js')

const addDAGNode = (taskState) => Object.assign({}, taskState)

const postValidation = (store) => (taskState) => new Promise((resolve, reject) => {
  store.dispatch(addOutput(taskState.uid, taskState))

  resolve({
    postValidation: ['addOutput']
  })
})

module.exports = postValidation
