'use strict'

const validateOutput = (taskState) => new Promise((resolve, reject) => {
  resolve({
    validations: {},
    validated: true
  })
})

module.exports = validateOutput
