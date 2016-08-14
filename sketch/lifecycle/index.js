'use strict'

const asyncDone = require('async-done')

// Lifecycle methods
const create = require('./create.js')
const resumable = require('./resumable.js')
const resolveInput = require('./resolve-input.js')
const createOperation = require('./create-operation.js')
const settleOperation = require('./settle-operation.js')
const resolveOutput = require('./resolve-output.js')
const validateOutput = require('./validate-output.js')
const postValidation = require('./post-validation.js')

const lifecycleMethods = {
  create,
  resumable,
  resolveInput,
  createOperation,
  settleOperation,
  resolveOutput,
  validateOutput,
  postValidation
}

module.exports = lifecycleMethods
