'use strict'

const { task } = require('.')

// const actions = require('./lib/get-actions.js')(store.dispatch)
//
// actions.setConfigValue('resume', 'on')
// actions.addTask({
//   name: 'foo'
// })
// actions.addTask({
//   name: 'bar',
//   params: { extreme: true }
// })

const intoStream = require('into-stream')
const fs = require('fs')

const myTask = task({
  input: { value: 'foo' },
  output: { file: 'file.txt' },
  name: 'My Task'
}, ({ input }) => intoStream(input).pipe(fs.createWriteStream('file.txt')))

myTask()
  .on('close' , function () {
    this.output()
  })
