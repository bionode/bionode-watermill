'use strict'

const Promise = require('bluebird')
const mkdirp = Promise.promisify(require('mkdirp'))

const setupTaskdir = (taskState, logger) => new Promise((resolve, reject) => {
  mkdirp(taskState.dir)
    .then(() => {
      logger.emit('log', `Made dir: ${taskState.dir}`)
      resolve(`Made dir: ${taskState.dir}`)
    })
    .catch((err) => {
      logger.emit('log', `Error creating ${taskState.dir}: ${err}`)
      return reject(`Error creating ${taskState.dir}: ${err}`)
    })
})

module.exports = setupTaskdir
