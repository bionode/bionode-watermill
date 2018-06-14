'use strict'

const { task, join } = require('../..')


const anotherTask = task({
  name: 'Run a JS file',
  output: 'ran.txt'
}, () => `node /home/evoxtorm/Desktop/Bionode-watermill/bionode-watermill/test/workflow-status/wait-error.js && touch ran.txt`)

anotherTask()



