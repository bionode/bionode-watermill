'use strict'

const globby = require('globby')
const fs = require('fs')

const { tab } = require('./utils.js')

// const existenceCheck = (obj) => {
//   const pattern = obj.file
//   const matches = globby.sync(pattern, obj.globOpts ? obj.globOpts : {})
//
//   if (matches.length === 0) {
//     obj.resolved = false
//
//     return false
//   } else {
//     obj.file = matches
//     obj.resolved = true
//
//     return true
//   }
// }

const isFileNull = (file) => {
  let stats
  try {
    stats = fs.statSync(file)
  } catch (err) {
    // TODO consider how else to handle this
    return false
    // console.log(err)
  }

  if (stats.size === 0) {
    return false
  } else {
    return true
  }
}

const nullCheck = (obj) => {
  if (!obj.resolved) {
    return false
  }

  obj.file.forEach((file) => {
    let notEmpty = isFileNull(file)

    if (!notEmpty) {
      return false
    }
  })

  return true
}

const matchTime = ({ file }) => {
  if (Array.isArray(file)) file = file[0]

  const fs = require('fs')

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
  isFileNull,
  nullCheck,
  matchTime
}
