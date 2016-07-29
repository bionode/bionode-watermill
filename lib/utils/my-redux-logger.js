'use strict'

const createLogger = require('redux-logger')
const chalk = require('chalk')

// Need to make a bunch of changes since this runs in regular terminal, not chrome
const logger = createLogger({
  colors: {
    title: false,
    prevState: false,
    action: false,
    nextState: false,
    error: false
  },
  logger: {
    log: (type, obj) => {
      let colour = 'white'
      if (/^action/.exec(type)) {
        if (/^action\ @/.exec(type)) {
          colour = 'red'
        } else {
          colour = 'blue'
        }
      } else if (/^prev/.exec(type)) {
        colour = 'grey'
      } else if (/^next/.exec(type)) {
        colour = 'magenta'
      }

      if (obj) {
        console.log(chalk[colour](type), JSON.stringify(obj, null, 2))
      } else {
        console.log(chalk[colour](type))
      }
    }
  },
  diff: true
})

module.exports = logger
