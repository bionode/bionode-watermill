'use strict'

// Lifecycle methods
const generateUid = require('./generate-uid.js')
const create = require('./create.js')
const checkResumable = require('./check-resumable.js')
const resolveInput = require('./resolve-input.js')
const createOperation = require('./create-operation.js')
const settleOperation = require('./settle-operation.js')
const resolveOutput = require('./resolve-output.js')
const validateOutput = require('./validate-output.js')
const postValidation = require('./post-validation.js')

const lifecycleMethods = {
  generateUid,
  create,
  checkResumable,
  resolveInput,
  createOperation,
  settleOperation,
  resolveOutput,
  validateOutput,
  postValidation
}

module.exports = lifecycleMethods
