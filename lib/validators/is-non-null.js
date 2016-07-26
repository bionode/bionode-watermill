'use strict'

const Promise = require('bluebird')
const fs = Promise.promisifyAll(require('fs'))

const isNonNull = (file) => new Promise((resolve, reject) => {
  return fs.statAsync(file)
  .then((stats) => {
    if (stats.size === 0) {
      reject(`File ${file} is null`)
    } else {
      resolve()
    }
  })
  .catch(reject)
})

module.exports = isNonNull
