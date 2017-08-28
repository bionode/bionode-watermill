'use strict'

const { assert } = require('chai')
const fs = require('fs')
const path = require('path')

const { task, join, junction, store } = require('../')

describe('junction', function() {
  it('join(junction(A, B), C) should work', function(done) {
    const simpleFileMaker = (content, ext) => task({
      output: `*.${ext}`,
      name: 'Simple file maker with ' + content
    }, ({ input, dir }) => {
      const ws = fs.createWriteStream(dir + '/' + content + '.' + ext)
      ws.write(content)
      ws.end()

      return ws
    })

    const concatFiles = task({
      input: {
        foo: '*.foo',
        bar: '*.bar'
      },
      name: 'Concat files',
      output: '*.foobar'
    }, ({ input, dir }) => {
      const ws = fs.createWriteStream(dir + '/' + 'concat.foobar')
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
      const lastTaskUid = results.context.trajectory[results.context.trajectory.length - 1]
      const finalOutput = store.getState().collection.vertexValue(lastTaskUid).resolvedOutput
      fs.readFile(finalOutput, (err, data) => {
        if (err) return done(err)
        assert.equal(data, 'foobar')
        done()
      })
    })

  })
})


