'use strict'

/**
 * Lifecycle method: create
 */

const { defaultTask, defaultTaskTypes } = require('../constants/default-task-state.js')
const stringify = require('json-stable-stringify')
const { hash } = require('../utils/index.js')

/**
 * Initialize a task object using defaultTask and user passed in options.
 *
 * @param props {Object}
 * @returns a task object
 */
const create = (props) => {
  // Check if props used a bad type for a certain key, if so, throw
  for (const key in defaultTaskTypes) {
    const val = defaultTaskTypes[key]

    const type = typeof(props[key])

    if (props[key] && type !== val) {
      throw new Error(`Bad type for ${key}, got ${type} but expected ${val}`)
    }
  }

  // TODO a "deepObjectAssign", so that properties like params can be extended

  const task = Object.assign({}, defaultTask, props, { created: Date.now() })
  const hashes = {
    input: hash(stringify(task.input)),
    output: hash(stringify(task.output)),
    params: hash(stringify(task.params))
  }
  const uid = hash(hashes.input + hashes.output + hashes.params)
  Object.assign(task, { hashes, uid })

  return task
}

module.exports = create
