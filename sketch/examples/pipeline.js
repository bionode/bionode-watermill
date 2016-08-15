'use strict'

const { task, join, parallel } = require('../')

const dumpPIDs = task({
  input: null,
  output: '*.pids',
  name: 'Dump all PIDs to *.pids'
}, () => `ps aux | awk '{print $2}' | tail -n +2 > ${Date.now()}.pids`)

dumpPIDs().then(console.log)
