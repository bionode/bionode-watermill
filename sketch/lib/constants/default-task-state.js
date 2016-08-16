'use strict'

const defaultConfig = require('./default-config-state.js')
const statusTypes = require('./task-status-types.js')

const defaultHashes = {
  input: null,
  output: null,
  params: null
}
exports.defaultHashes = defaultHashes

exports.defaultTask = {
  threads: defaultConfig.threads,
  container: defaultConfig.container,
  resume: defaultConfig.resume,
  uid: null,
  hashes: defaultHashes,
  name: 'Unnamed Task',
  dir: null,
  input: null,
  output: null,
  operation: null,
  status: statusTypes.STATUS_CREATED,
  created: null,
  validated: false,
  params: {},
  trajectory: []
}

exports.defaultTaskTypes = {
  threads: 'number',
  container: 'string',
  resume: 'string', // TODO 'on' or 'off'
  uid: 'string',
  hashes: 'object',
  name: 'string',
  // TODO rest
}

