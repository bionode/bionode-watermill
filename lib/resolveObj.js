'use strict'

const { tab } = require('./utils.js')
const chalk = require('chalk')


// TODO error if value, file, etc, not at bottom of object
const resolveObj = (obj) => {
  if (obj instanceof Array) {
    return obj.map(item => resolveObj(item))
  }

  for (let key in obj) {
    let val = obj[key]

    switch (key) {
      case 'value':
        // Simply resolve { value: 'foo' } to 'foo'
        obj = val 
        console.log(tab(2) + 'input: ' + chalk.magenta(val) + ' from ' + chalk.cyan('value'))
        break
      default:
        console.log(val)
        console.log(obj)
        console.log('oops')
        // obj[key] = resolveObj(val)
    } 
  }

  return obj
}

module.exports = resolveObj
