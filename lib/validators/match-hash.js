'use strict'

const fs = require('fs')
const checksum = require('checksum')
const { tab } = require('../utils')

const matchHash = (file) => new Promise((resolve, reject) => {
  let currentLogs
  try {
    currentLogs = JSON.parse(fs.readFileSync('watermill.json'))
  } catch (err) {
    reject('No watermill.json file')
  }

  if (!currentLogs[file]) {
    reject('File not in watermill.json')
  }

  checksum.file(file, (err, sum) => {
    if (err) {
      reject(err)
    }

    if (!currentLogs[file]) {
      return reject(`No key for ${file} in watermill.json`)
    }

    if (currentLogs[file].hash === sum) {
      resolve(`Hash matched for ${file}`)
    } else {
      reject(`Hash did not match for ${file}`)
    }
  })
})

module.exports = matchHash
