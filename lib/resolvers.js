'use strict'

// pattern matching
const globby = require('globby')
const multimatch = require('multimatch')
const chalk = require('chalk')

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

const fileResolver = (key, obj, mode, input) => {
  const val = obj[key]

  if (mode === 'passed') {
    // Compare to task input
    // Can also check obj.resolved here
    // TODO handle better input object 
    const matches = multimatch([val], [input.file])

    if (matches.length === 1) {
      console.log(tab(2) + 'input: ' + chalk.magenta(val) + ' from passed in ' + chalk.cyan('file'))
      obj = val
    }
  } else if (mode === 'start') {
    const matches = globby.sync(val)

    if (matches.length >= 0) {
      console.log(tab(2) + 'input: ' + chalk.magenta(matches[0]) + ' fromm ' + chalk.cyan('file'))
      obj = matches[0]
    }
  } else if (mode === 'end') {
    console.log(tab(2) + chalk.yellow('Looking') + ' for ' + val)
    const globbyData = globby.sync(val)

    let newObj
    if (globbyData.length > 0) {
      console.log(tab(3) + chalk.green(' found ') + globbyData[0])
      newObj = { file: globbyData[0], resolved: true }
    } else {
      newObj = Object.assign(val, { resolved: false })
    }

    obj = newObj

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
