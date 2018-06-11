'use strict'

const { task, join } = require('../../..')

const fs = require('fs')
const through = require('through2')
const { execFile } = require('child_process')

// const myTask = task({
//   input: '*.source.txt',
//   output: '*.sink.txt',
//   name: 'Capitalize alphabet'
// }, ({ input }) => {
//     fs.createReadStream(input)
//       .pipe(through(function (chunk, enc, next) {
//         this.push(chunk.toString().toUpperCase())

//         next()
//       }))
//       .pipe(fs.createWriteStream(input.split('/').slice(0, -1).join('/') + '/alphabet.sink.txt'))
//   }
// )

const anotherTask = task({
  name: 'Run a JS file',
  output: 'ran.txt'
}, () => `node wait-error.js && touch ran.txt`)


anotherTask()

// const pipeline = join(myTask, anotherTask)
// pipeline()


