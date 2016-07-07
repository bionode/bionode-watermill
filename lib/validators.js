'use strict'

const globby = require('globby')
const fs = require('fs')


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

const nullCheck = (file) => {
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

module.exports = {
  existenceCheck,
  nullCheck
}
