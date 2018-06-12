'use strict'

const { task, join } = require('../../..')

const fs = require('fs')
const through = require('through2')
const { execFile } = require('child_process')


const anotherTask = task({
  name: 'Run a JS file',
  output: 'ran.txt'
}, () => `node wait-error.js && touch ran.txt`)


anotherTask()



