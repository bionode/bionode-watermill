'use strict'

const crypto = require('crypto')

const tab = (num, size = 2) => {
  let tab = ''
  for (let i=0; i < size; i++) {
    tab += ' '
  }

  let str = ''
  for (let i=0; i < num; i++) {
    str += tab
  }

  return str
}

const isChildProcess = (stream) => stream.stdio && Array.isArray(stream.stdio) && stream.stdio.length === 3

const hash = (data) => crypto.createHash('sha256').update(data).digest('hex')

module.exports = {
  tab,
  isChildProcess,
  hash
}
