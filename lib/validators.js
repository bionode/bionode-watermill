'use strict'

const globby = require('globby')

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

module.exports = {
  existenceCheck
}
