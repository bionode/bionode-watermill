'use strict'

const addDAGNode = (taskState) => Object.assign({}, taskState)

const postValidation = (taskState) => new Promise((resolve, reject) => {
  resolve({
    postValidation: {}
  })
})

module.exports = postValidation
