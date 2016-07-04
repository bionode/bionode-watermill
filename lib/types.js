'use strict'

// mostly for instanceof checks inside task
// TODO refactor out for vanilla objs
class File {
  constructor(value) {
    this.value = value
    this.resolved = null

    this.resolve = this.resolve.bind(this)
  }

  resolve(value) {
    this.resolved = value
  }
}

exports.File = File

const types = [
  'file',
  'value',
  'stream'
]

exports.types = types
