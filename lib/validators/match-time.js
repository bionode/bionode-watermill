'use strict'

const globby = require('globby')
const fs = require('fs')

const { tab } = require('../utils.js')

const matchTime = (file) => {
  let currentLogs
  try {
    let content = fs.readFileSync('waterwheel.json')
    currentLogs = JSON.parse(content)
  } catch(err) {
    console.log('errored')
    return false
  }

  if (!currentLogs[file]) {
    return false
  }

  let stats
  try {
    stats = fs.statSync(file)
  } catch(err) {
    return false
  }

  const time = stats.mtime.getTime()

  if (currentLogs[file].time === time) {
    console.log(tab(3) + 'Time matched for ' + file)
    return true
  } else {
    console.log(tab(3) + 'Time did not match for ' + file)
    return false
  }
}

module.exports = {
  matchTime
}
