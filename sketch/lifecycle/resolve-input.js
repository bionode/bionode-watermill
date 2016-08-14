'use strict'

/**
 * Lifecycle method: resolve input
 */

const resolveInput = (taskState) => {
  return {
    resolvedInput: taskState.input
  }
}

module.exports = resolveInput
