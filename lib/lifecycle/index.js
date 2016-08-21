'use strict'

// Lifecycle methods
const generateUid = require('./generate-uid.js')
const checkResumable = require('./check-resumable.js')
const resolveInput = require('./resolve-input.js')
const setupTaskdir = require('./setup-taskdir.js')
const createOperation = require('./create-operation.js')
const settleOperation = require('./settle-operation.js')
const resolveOutput = require('./resolve-output.js')
const validateOutput = require('./validate-output.js')
const postValidation = require('./post-validation.js')

const lifecycleMethods = {
  generateUid,
  checkResumable,
  resolveInput,
  setupTaskdir,
  createOperation,
  settleOperation,
  resolveOutput,
  validateOutput,
  postValidation
}

module.exports = lifecycleMethods
