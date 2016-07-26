'use strict'

const { task, shell } = require('../..')

const intoStream = require('into-stream')
const fs = require('fs')
const through = require('through2')

const throughUpperCase  = through(function (chunk, env, next) {
  this.push(chunk.toString().toUpperCase())

  next()
})

const myTask = task({
  input: 'sources/*.source.txt',
  output: '*.sink.txt',
  name: 'Capitalize alphabet'
}, ({ input }) => {

  let stream
  for (let item of input) {
    const name = item.split('/').pop().split('.')[0]

    stream = fs.createReadStream(item)
      .pipe(throughUpperCase)
      .pipe(fs.createWriteStream(name + '.sink.txt'))
  }

  // return last stream...
  return stream
})

myTask()
  .on('close' , function () {
    console.log('============== GOT CLOSE =======================')
    // creating output on close
    // this.output()

    // using already created output (preferred, but child process throws close itself)
    console.log('output:', this._output)
  })
