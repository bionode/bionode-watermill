'use strict'

// pattern matching
const globby = require('globby')
const multimatch = require('multimatch')
const chalk = require('chalk')

const validators = require('./validators.js')
const { tab } = require('./utils.js')
const traverser = require('./traverser.js')

/**
 * key = value
 * obj = { value: 'foo' }
 * @returns 'foo'
 */
const valueResolver = (key, obj, mode) => {
  const val = obj[key]

  console.log(tab(2) + 'resolved ' + chalk.magenta(val) + (mode === 'passed' ? ' from passed in ' : ' from ') + chalk.cyan('value'))
  return val
}

const defaultFileValidators = [
  validators.existenceCheck,
  validators.nullCheck
]

const fileResolver = (key, obj, mode, input) => {
  const val = obj[key]

  if (mode === 'passed') {
    // Compare to task input
    // Can also check obj.resolved here
    // TODO handle better input object
    const matches = multimatch(val, [input.file])

    // Just check if patterns match. Validation occurred after end of last task.

    if (matches.length === 1) {
      console.log(tab(2) + 'input: ' + chalk.magenta(val) + ' from passed in ' + chalk.cyan('file'))
      obj = val[0]
    }
  } else if (mode === 'start') {
    defaultFileValidators.forEach((validator) => {
      const result = validator(obj)

      if (!result) {
        console.log(tab(2) + chalk.red('An input validation failed'))
        return obj
      }
    })

    // If we got here, no validators returned false

    console.log(tab(2) + 'input: ' + chalk.magenta(obj.file[0]) + ' from ' + chalk.cyan('file'))
    obj = obj.file[0]
  } else if (mode === 'end') {
    console.log(tab(2) + chalk.yellow('Running validations') + ' for ' + val)
    defaultFileValidators.forEach((validator) => {
      const result = validator(obj)

      if (!result) {
        console.log(tab(2) + chalk.red('An input validation failed'))
        return obj
      }
    })
    console.log(tab(3) + '👌  passed validators')
    // TODO check for resolved:true
  }

  return obj
}

const outerResolver = (obj, props, mode) => {
  return traverser(obj, (key, obj) => {
      const val = obj[key]

      switch(key) {
        case 'value':
          obj = valueResolver(key, obj, mode)
          break
        case 'file':
          obj = fileResolver(key, obj, mode, props.input)
          break
        case 'resolved':
          break
        default:
          console.log(key, val)
          return false
      }

      return obj
    }) 
}

module.exports = outerResolver
