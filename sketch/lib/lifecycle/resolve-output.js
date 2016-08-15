'use strict'

const _ = require('lodash')

const matchToFs = require('../matchers/match-to-fs.js')
const applicator = require('../utils/applicator.js')

const resolveOutput = (taskState) => new Promise((resolve, reject) => {
  const { output, dir } = taskState

  let opts = {}
  if (!_.isNull(dir)) {
    opts.cwd = dir
  }

  applicator(output, (item) => matchToFs(item, null, opts))
  .then(results => {
    resolve({ success: true, resolvedOutput: results })
  })
  // Output was not resolved
  .catch((err) => {
    reject(err)
  })

})

module.exports = resolveOutput
