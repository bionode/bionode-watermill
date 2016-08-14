'use strict'

const resolveOutput = (taskState) => new Promise((resolve, reject) => {
  resolve({ resolvedOutput: taskState.output })
})

module.exports = resolveOutput
