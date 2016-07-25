'use strict'

const { task, shell } = require('.')

const intoStream = require('into-stream')
const fs = require('fs')
const myTask = task({
  input: { value: 'foo' },
  output: '*ile.txt',
  name: 'My Task'
// }, ({ input }) => shell(`echo "${input}\n${input}" > file.txt`))
}, ({ input }) => intoStream(input).pipe(fs.createWriteStream('file.txt')))

myTask()
  .on('close' , function () {
    console.log('============== GOT CLOSE =======================')
    // creating output on close
    // this.output()

    // using already created output (preferred, but child process throws close itself)
    console.log('output:', this._output)
  })
