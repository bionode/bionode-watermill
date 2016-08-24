'use strict'

const globby = require('globby')
const fs = require('fs')

const matchTime = (file) => new Promise((resolve, reject) => {
  let currentLogs
  try {
    let content = fs.readFileSync('watermill.json')
    currentLogs = JSON.parse(content)
  } catch(err) {
    console.log('errored')
    reject()
  }

  if (!currentLogs[file]) {
    reject('Path not in watermill.json')
  }

  let stats
  try {
    stats = fs.statSync(file)
  } catch(err) {
    reject()
  }

  const time = stats.mtime.getTime()

  if (currentLogs[file].time === time) {
    console.log('Time matched for ' + file)
    resolve()
  } else {
    reject('Time did not match for ' + file)
  }
})

module.exports = matchTime
