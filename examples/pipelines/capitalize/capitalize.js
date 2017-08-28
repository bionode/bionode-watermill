'use strict'

const { task } = require('../../..')

const fs = require('fs')
const through = require('through2')

const myTask = task({
  input: '*.source.txt',
  output: '*.sink.txt',
  name: 'Capitalize alphabet'
}, ({ input }) => {
    fs.createReadStream(input)
      .pipe(through(function (chunk, enc, next) {
        this.push(chunk.toString().toUpperCase())

        next()
      }))
      .pipe(fs.createWriteStream(input.split('/').slice(0, -1).join('/') + '/alphabet.sink.txt'))
  }
)

myTask()
