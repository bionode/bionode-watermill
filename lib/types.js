'use strict'

// mostly for instanceof checks inside task
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
