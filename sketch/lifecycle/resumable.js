'use strict'

const resumable = (taskState) => {
  return taskState.resumable === 'on' ? true : false
}

module.exports = resumable
