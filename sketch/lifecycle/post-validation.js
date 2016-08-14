'use strict'

const addDAGNode = (taskState) => Object.assign({}, taskState)

const postValidation = (taskState) => {
  return addDAGNode(taskState)
}

module.exports = postValidation
