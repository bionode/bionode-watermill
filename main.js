'use strict'

const { task, shell } = require('.')

const intoStream = require('into-stream')
const fs = require('fs')
const through = require('through2')

const myTask = task({
  input: '*.source.txt',
  output: '*.sink.txt',
  name: 'Capitalize alphabet'
}, ({ input }) =>
  fs.createReadStream(input)
    .pipe(through(function (chunk, enc, next) {
      this.push(chunk.toString().toUpperCase())

      next()
    }))
    .pipe(fs.createWriteStream('alphabet.sink.txt'))
)

myTask()
  .on('close' , function () {
    console.log('============== GOT CLOSE =======================')
    // creating output on close
    // this.output()

    // using already created output (preferred, but child process throws close itself)
    console.log('output:', this._output)
  })
