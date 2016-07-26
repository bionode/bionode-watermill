'use strict'

const globby = require('globby')
const fs = require('fs')

const { tab } = require('../utils.js')

const matchTime = (file) => new Promise((resolve, reject) => {
  let currentLogs
  try {
    let content = fs.readFileSync('waterwheel.json')
    currentLogs = JSON.parse(content)
  } catch(err) {
    console.log('errored')
    reject()
  }

  if (!currentLogs[file]) {
    reject('Path not in waterwheel.json')
  }

  let stats
  try {
    stats = fs.statSync(file)
  } catch(err) {
    reject()
  }

  const time = stats.mtime.getTime()

  if (currentLogs[file].time === time) {
    console.log(tab(3) + 'Time matched for ' + file)
    resolve()
  } else {
    reject(tab(3) + 'Time did not match for ' + file)
  }
})

module.exports = matchTime
