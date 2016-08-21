'use strict'

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

module.exports = {
  tab,
  isChildProcess
}
