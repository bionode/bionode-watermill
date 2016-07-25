'use strict'

const Promise = require('bluebird')
const fs = Promise.promisifyAll(require('fs'))
const path = require('path')
const { assert } = require('chai')

const twd = path.resolve(__dirname, 'files', 'output-resolver')

const globby = require('globby')
const chalk = require('chalk')
const _ = require('lodash')

const matchToFs = (pattern, num = 1, opts = {}) => new Promise((resolve, reject) => {
  if (_.isNull(num)) num = 1

  let _num
  if (_.isNumber(num)) _num = (len) => len === num
  if (_.isFunction(num)) _num = num

  globby(pattern, opts)
  .then((paths) => {
    // You could expect 0 matches with a custom `num` function if desired
    if (_num(paths.length)) {
      // Turn one item arrays into just the item
      if (paths.length === 1) paths = paths[0]
      // Else, resolve with array of paths
      return resolve(paths)
    } else {
      // TODO more descriptive error message if num was a number
      reject(`No matches found for pattern ${chalk.magenta(pattern)} in ${opts.cwd || process.cwd()}`)
    }
  })
  .catch(reject)
})


const resolver = (output, matcher) => {
  if (_.isString(output)) {
    return matcher(output)
  }
}

describe('Output resolver', function() {
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
})
