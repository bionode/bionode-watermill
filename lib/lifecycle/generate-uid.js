'use strict'

const stringify = require('json-stable-stringify')

const hash = require('../utils/hash.js')

const { defaultHashes } = require('../constants/default-task-state.js')

/**
 * Generate a task UID and hashes for input, output, params
 *
 * @param {Object} desctructure out input, output, params
 *
 * @returns {Object} of the structure { uid, hashes: { input, output, params } }
 */
const generateUid = ({ input, output, params, name }) => {
  const userIOP = { input, output, params, name }

  // Clean out undefined or null values
  // Thus is the user passes input, output, but no params
  // params is still set to the default instead of undefined
  const cleanIOP = {}
  for (const key in userIOP) {
    if (userIOP[key]) cleanIOP[key] = userIOP[key]
  }

  // Assign user passed in values over defaults
  const IOP = Object.assign({}, defaultHashes, cleanIOP)

  // Map each value in IOP to its hash
  // Use json-stable-stringify which is deterministic
  const hashes = Object.keys(IOP).reduce(
    (sum, key) => Object.assign({}, sum, { [key]:hash(stringify(IOP[key])) }),
    {}
  )

  // uid is a hash of hashes
  const uid = hash(hashes.input + hashes.output + hashes.params + hashes.name)

  const result = {
    uid,
    hashes
  }

  return result
}

module.exports = generateUid
