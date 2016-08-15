'use strict'

/**
 * Lifecycle method: resolve input
 */

const resolveInput = (taskState) => new Promise((resolve, reject) => {
  resolve({
    resolvedInput: taskState.input
  })
})

module.exports = resolveInput
