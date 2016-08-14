'use strict'

/**
 * Lifecycle method: create
 */

/**
 * Create first taskState given props
 */
const create = (props) => {
  return Object.assign({}, props, { created: Date.now() })
}

module.exports = create
