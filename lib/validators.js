'use strict'

const globby = require('globby')
const fs = require('fs')


const existenceCheck = (obj) => {
  const pattern = obj.file
  const matches = globby.sync(pattern, obj.globOpts ? obj.globOpts : {})

  if (matches.length === 0) {
    obj = Object.assign(obj, { resolved: false })
  } else {
    obj = Object.assign(obj, { file: matches }, { resolved: true })
  }

  return obj
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
