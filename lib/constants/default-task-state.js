'use strict'

const defaultConfig = require('./default-config-state.js')
const statusTypes = require('./task-status-types.js')

const defaultHashes = {
  input: null,
  output: null,
  params: null
}
exports.defaultHashes = defaultHashes

const defaultContext = {
  trajectory: []
}
exports.defaultContext = defaultContext

const defaultLog = {
  log: Buffer.from(''),
  // Some of these we dont want to store - e.g. some programs use stderr for
  // logging, stdout for the real data, (and dont always exit with proper error
  // codes..)
  stdout: Buffer.from(''),
  stderr: Buffer.from(''),
  currentLog: Buffer.from(''),
  // will need to copy into a new buffer each action as per state immutability
  exitCode: null,
}

exports.defaultTask = {
  threads: defaultConfig.threads,
  container: defaultConfig.container,
  resume: defaultConfig.resume,
  uid: null,
  instance: null,
  hashes: defaultHashes,
  name: 'Unnamed Task',
  dir: defaultConfig.workdir,
  input: null,
  output: null,
  operation: null,
  operationString: null,  // just one is really necessary between this and
  // task.js definition
  status: statusTypes.STATUS_CREATED,
  created: null,
  validated: false,
  params: {},
  // trajectory: [],
  context: defaultContext,
  log: defaultLog,
  type: 'task'
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

