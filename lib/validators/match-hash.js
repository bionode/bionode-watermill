'use strict'

const fs = require('fs')
const checksum = require('checksum')
const { tab } = require('../utils')

const matchHash = (file) => new Promise((resolve, reject) => {
  let currentLogs
  try {
    currentLogs = JSON.parse(fs.readFileSync('waterwheel.json'))
  } catch (err) {
    reject('No waterwheel.json file')
  }

  if (!currentLogs[file]) {
    reject('File not in waterwheel.json')
  }

  checksum.file(file, (err, sum) => {
    if (err) {
      reject(err)
    }

    if (currentLogs[file].hash === sum) {
      resolve(`Hash matched for ${file}`)
    } else {
      reject(`Hash did not match for ${file}`)
    }
  })
})

module.exports = matchHash
