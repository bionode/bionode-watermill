'use strict'

const { task, join } = require('../..')
const { waiterror } = require('./wait-error.js')


const anotherTask = task({
  name: 'Run a JS file',
  output: 'ran.txt'
}, () => `node wait-error.js && touch ran.txt`)

anotherTask()



