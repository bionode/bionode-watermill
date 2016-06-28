'use strict'

const stream = require('readable-stream')

class Duplex extends stream.Duplex {
  constructor(opts) {
    super(opts)
  }

  _read(size) {
    console.log('reading')
  }

  _write(chunk, encoding, callback) {
    console.log('writing: ' + chunk.toString())

    this.push(chunk)
    callback()
  }
}

// Would use https://www.npmjs.com/package/stream-from-promise but may
// need to handle special resolve values (e.g. promise resolves to a stream) later
class ReadableFromPromise extends stream.Readable {
  constructor(promise, opts) {
    super(opts)

    this._promise = promise
  }

  _read(size) {
    this._promise
      .then((data) => {
        this.push(data)
        this.push(null)
      })
      .catch((err) => {
        this.emit('error', err)
        this.push(null)
      })
  }
}

/**
 * The Task class.
 * @param  {Object} props         object with input and output
 * @param  {function} actionCreator function to produce an actionCreator
 * @return {stream}               stream to be orchestrated
 */
const task = (props, actionCreator) => {
  // TODO check if actionCreator is a function

  // TODO
  // 1. resolve props OR
  // 2. use values passed in from last Task

  // TODO use input/output, resolved promise, callback to determine stream type.
  // Or just always be duplex.
  let stream
  // const stream = new Readable()

  // Create the action
  const action = actionCreator(props)

  // Get its type: can be stream, promise, function
  const actionType = (() => {
    if (action instanceof Promise) {
      return 'promise'
    } else {
      return 'unkown'
    }
  })()

  switch(actionType) {
    case 'promise':
      stream = new ReadableFromPromise(action)
      break
    default:
      console.log('Bad action type')
  }

  return stream
}

module.exports = task
