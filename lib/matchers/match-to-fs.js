'use strict'

const _ = require('lodash')
const globby = require('globby')
const chalk = require('chalk')

const matchToFs = (pattern, num = 1, opts = {}) => new Promise((resolve, reject) => {
  if (_.isNull(num) || _.isUndefined(num)) num = (len) => len !== 0

  let _num
  if (_.isNumber(num)) _num = (len) => len === num
  if (_.isFunction(num)) _num = num

  opts = Object.assign({}, { realpath: true }, opts)

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
  .catch((err) => {
    console.log('math-to-fs error: ', err)
    reject(err)
  })
})

module.exports = matchToFs
