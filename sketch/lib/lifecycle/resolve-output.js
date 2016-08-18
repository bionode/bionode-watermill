'use strict'

const _ = require('lodash')

const matchToFs = require('../matchers/match-to-fs.js')
const applicator = require('../utils/applicator.js')

const resolveOutput = (taskState, logger) => new Promise((resolve, reject) => {
  const { uid, output, dir } = taskState
  logger.emit('log', `Starting resolve output for ${uid}`)

  let opts = {}
  if (!_.isNull(dir)) {
    opts.cwd = dir
  }

  applicator(output, (item) => matchToFs(item, null, opts))
    .then(results => resolve(results))
    .catch(err => reject(err))
})

module.exports = resolveOutput
