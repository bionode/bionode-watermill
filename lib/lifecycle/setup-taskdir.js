'use strict'

const Promise = require('bluebird')
const mkdirp = Promise.promisify(require('mkdirp'))
const path = require('path')
const link = require('fs-symlink')

const applicator = require('../utils/applicator.js')

const setupTaskdir = (taskState, logger) => new Promise((resolve, reject) => {
  // Make the task directory
  mkdirp(taskState.dir)
    .then(() => {
      logger.emit('log', `Made dir: ${taskState.dir}`)
      // resolve(`Made dir: ${taskState.dir}`)
    })
    .catch((err) => {
      logger.emit('log', `Error creating ${taskState.dir}: ${err}`)
      return reject(`Error creating ${taskState.dir}: ${err}`)
    })

  // Create symlinks for resolved inputs
  const symlinkToDir = (source) => new Promise((resolve, reject) => {
    console.log('Gonna make a symlink for: ', source)
    const dest = path.resolve(taskState.dir, source.split('/').pop())
    console.log('Into: ', dest)

    link(source, dest, 'junction').then(() => {
      resolve(dest)
    }).catch(err => console.log('symlink error: ', err))
  })

  applicator(taskState.resolvedInput, symlinkToDir).then((results) => {
    console.log('symlinkedInput: ', results)
    resolve(results)
  })
})

module.exports = setupTaskdir
