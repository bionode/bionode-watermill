'use strict'

const { task, shell } = require('../../..')()

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
  .on('task.finish' , function (task) {
    console.log('output:', task.resolvedOutput)
  })
