'use strict'

const stream = require('readable-stream')

class StreamFromCallback extends stream.Readable {
  constructor(action, opts) {
    super(opts)

    this._action = action
  }

  _read(size) {
    this._action((err, data) => {
      if (err !== null) {
        this.emit('error', err)
        this.push(null)
      } else {
        this.push(data)
        this.push(null)
      }
    })
  }
}

module.exports = StreamFromCallback
