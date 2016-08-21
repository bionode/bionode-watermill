'use strict'

const { assert } = require('chai')
const fs = require('fs')
const path = require('path')

const twd = path.resolve(__dirname, 'files', 'junction')

const { task, join, junction } = require('../')

describe('junction', function() {
  it('join(junction(A, B), C) should work', function(done) {
    const simpleFileMaker = (content, ext) => task({
      output: `${twd}/*.${ext}`,
      name: 'Simple file maker with ' + content
    }, ({ input }) => {
      const ws = fs.createWriteStream(twd + '/' + content + '.' + ext)
      ws.write(content)
      ws.end()

      return ws
    })

    const concatFiles = task({
      input: {
        foo: `${twd}/*.foo`,
        bar: `${twd}/*.bar`
      },
      name: 'Concat files',
      output: '${twd}/*.foobar'
    }, ({ input }) => {
      const ws = fs.createWriteStream(twd + '/' + 'concat.foobar')
      const rsFoo = fs.createReadStream(input.foo)
      const rsBar = fs.createReadStream(input.bar)

      rsFoo.pipe(ws, { end: false })
      rsFoo.on('end', () => rsBar.pipe(ws))

      return ws
    })

    const fooTask = simpleFileMaker('foo', 'foo')
    const barTask = simpleFileMaker('bar', 'bar')

    /* Looks like:
     * fooTask -->|
     *            |--> concatFiles
     * barTask -->|
     */
    const pipeline = join(
      junction(fooTask, barTask),
      concatFiles
    )

    pipeline().then((results) => {
      console.log('pipeline: ', results)
      fs.readFile(twd + '/' + 'concat.foobar', (err, data) => {
        if (err) return done(err)
        assert.equal(data, 'foobar')
        done()
      })
    })

  })
})


