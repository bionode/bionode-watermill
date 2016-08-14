'use strict'

const addDAGNode = (taskState) => Object.assign({}, taskState)

const postValidation = (taskState) => {
  return {
    postValidation: {}
  }
}

module.exports = postValidation
