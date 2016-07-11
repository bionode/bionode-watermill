'use strict'

const globby = require('globby')
const fs = require('fs')

const { tab } = require('./utils.js')


const existenceCheck = (obj) => {
  const pattern = obj.file
  const matches = globby.sync(pattern, obj.globOpts ? obj.globOpts : {})

  if (matches.length === 0) {
    obj.resolved = false
    
    return false
  } else {
    obj.file = matches
    obj.resolved = true
    
    return true
  }
}

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

const matchHash = ({ file }) => {
  if (Array.isArray(file)) file = file[0]

  const fs = require('fs')
  const checksum = require('checksum')

  let currentLogs
  try {
    currentLogs = JSON.parse(fs.readFileSync('waterwheel.json'))
  } catch (err) {
    return false
  }

  if (!currentLogs[file]) {
    return false
  }

  checksum.file(file, (err, sum) => {
    if (err) {
      return false
    }

    console.log((currentLogs[file].hash === sum ? 'hash matched' : 'hash didnt match') + ' for ' + file) 
  })

  // TODO async validators
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
  existenceCheck,
  isFileNull,
  nullCheck,
  matchHash,
  matchTime
}
