'use strict'

const Promise = require('bluebird')
const fs = Promise.promisifyAll(require('fs'))
const path = require('path')
const { assert } = require('chai')

const twd = path.resolve(__dirname, 'files', 'output-resolver')

const matchToFs = require('../lib/matchers/match-to-fs.js')
const resolver = require('../lib/resolver2.js')

describe('Resolver', function() {
  it('should resolve a single, hardcoded file', function(done) {
    const output = 'foo.txt'

    resolver(output, (item) => matchToFs(item, null, { cwd: twd }))
    .then((resolvedOutput) => {
      assert(resolvedOutput === 'foo.txt')

      done()
    })
    .catch(done)
  })

  it('should resolve a single, patterned file', function(done) {
    const output = '*_suffix.txt'

    resolver(output, (item) => matchToFs(item, 1, { cwd: twd }))
    .then((resolvedOutput) => {
      assert(resolvedOutput === 'something_suffix.txt')

      done()
    })
    .catch(done)
  })

  it('should reject if there is no match', function(done) {
    const output = 'nonexistent.txt'

    resolver(output, (item) => matchToFs(item, 1, { cwd: twd }))
    .then(() => done(new Error('Should not have resolved')))
    .catch((err) => done())
  })

  it('should reject if there are an unexpected amount of matches', function(done) {
    const output = '*_two.txt'

    resolver(output, (item) => matchToFs(item, 1, { cwd: twd }))
    .then(() => done(new Error('Should not have resolved')))
    .catch(err => done())
  })

  it('should resolve if there is an expected number of matches', function(done) {
    const output = '*_two.txt'

    resolver(output, (item) => matchToFs(item, (len) => len === 2, { cwd: twd }))
    .then((resolvedOutput) => {
      assert.deepEqual(resolvedOutput, ['a_two.txt', 'b_two.txt'], 'Did not produce expected matches')

      done()
    })
    .catch(done)
  })

  it('should resolve an array of patterns', function(done) {
    const output = ['*_suffix.txt', '*_two.txt']

    resolver(output, (item) => matchToFs(item, (len) => len !== 0, { cwd: twd }))
    .then((resolvedOutput) => {
      assert.deepEqual(resolvedOutput,
                       [ 'something_suffix.txt', ['a_two.txt', 'b_two.txt'] ],
                       'Did not produce expected matches')

      done()
    })
    .catch(done)
  })

  it('should reject if any of the patterns provided do not match', function(done) {
    const output = ['*_suffix.txt', 'nonexistent.txt']

    resolver(output, (item) => matchToFs(item, (len) => len !== 0, { cwd: twd }))
    .then(() => done(new Error('Should not have resolved')))
    .catch(err => done())
  })

  it('should handle simple objects', function(done) {
    const output = { a: 'a*.txt', b: 'b*.txt' }

    resolver(output, (item) => matchToFs(item, (len) => len !== 0, { cwd: twd }))
    .then((resolvedOutput) => {
      assert.deepEqual(resolvedOutput, {
        a: 'a_two.txt',
        b: 'b_two.txt'
      }, 'Resolved output not as expected')

      done()
    })
    .catch(done)
  })

  it('should handle arrays made of objects', function(done) {
    const output = [{suffixes: '*_suffix.txt'}, { twos: '*_two.txt' }]

    resolver(output, (item) => matchToFs(item, null, { cwd: twd }))
    .then((resolvedOutput) => {
      assert.deepEqual(resolvedOutput, [
        { suffixes: 'something_suffix.txt' },
        { twos: [ 'a_two.txt', 'b_two.txt' ] }
      ], 'Resolved output not as expected')

      done()
    })
    .catch(done)
  })

  it('should handle objects with values as arrays', function(done) {
    const output = { files: ['*_suffix.txt', '*.txt'], foo: 'foo.txt' }

    resolver(output, (item) => matchToFs(item, null, { cwd: twd }))
    .then((resolvedOutput) => {
      assert.deepEqual(resolvedOutput, {
        files: ['something_suffix.txt', ['a_two.txt', 'b_two.txt', 'foo.txt', 'something_suffix.txt']],
        foo: 'foo.txt'
      })

      done()
    })
    .catch(done)
  })

  it('should handle a complex mix of nested objects and arrays', function(done) {
    const output = {
      suffixes: {
        patterns: [{
          broad: [
            '*.txt',
            '*_two.txt'
          ]
        }, {
          fixed: '*_suffix.txt'
        }]
      }
    }

    resolver(output, (item) => matchToFs(item, null, { cwd: twd }))
    .then((resolvedOutput) => {
      assert.deepEqual(resolvedOutput, {
        suffixes: {
          patterns: [{
            broad: [
              ['a_two.txt', 'b_two.txt', 'foo.txt', 'something_suffix.txt'],
              ['a_two.txt', 'b_two.txt']
            ]
          }, {
            fixed: 'something_suffix.txt'
          }]
        }
      })

      done()
    })
    .catch(done)
  })
})
