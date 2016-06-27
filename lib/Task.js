'use strict'

const stream = require('readable-stream')

class Duplex extends stream.Duplex {
  constructor(opts) {
    super(opts)
  }

  _read() {
    console.log('reading?')
  }

  _write(chunk, encoding, callback) {
    console.log('writing: ' + chunk.toString())

    this.push(chunk)
    callback()
  }
}

/**
 * The Task class.
 * @param  {Object} props         object with input and output
 * @param  {function} actionCreator function to produce an actionCreator
 * @return {stream}               stream to be orchestrated
 */
const task = (props, actionCreator) => () => {
  // TODO check if actionCreator is a function

  // TODO
  // 1. resolve props OR
  // 2. use values passed in from last Task

  // TODO use input/output, resolved promise, callback to determine stream type.
  // Or just always be duplex.
  const stream = new Duplex()

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
      action
        .then(data => stream.end(data))
        .catch(err => console.error(err))

      break
    default:
      console.log('Bad action type')
  }

  return stream
}

module.exports = task
