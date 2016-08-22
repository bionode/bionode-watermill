'use strict'

const _ = require('lodash')

const matchToFs = require('../matchers/match-to-fs.js')
const applicator = require('../utils/applicator.js')

const resolveOutput = (taskState, logger) => new Promise((resolve, reject) => {
  // const { uid, output, dir } = taskState
  const { uid, output } = taskState
  const dir = taskState.dir
  logger.emit('log', `Starting resolve output for ${uid}`)

  if (!output) {
    logger.emit('log', 'Skipping since output is null', 1)
    return resolve()
  }

  let opts = {}
  if (!_.isNull(dir)) {
    opts.cwd = dir
  }

  applicator(output, (item) => matchToFs(item, null, opts))
    .then(results => resolve(results))
    .catch(err => reject(err))
})

module.exports = resolveOutput
