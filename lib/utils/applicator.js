'use strict'

const Promise = require('bluebird')
const _ = require('lodash')

const applicator = (obj, matcher) => {
  if (_.isArray(obj)) {
    return Promise.all(obj.map(item => applicator(item, matcher)))
  }

  if (_.isPlainObject(obj)) {
    return Promise.reduce(Object.keys(obj), (sum, currentKey) => {
      return applicator(obj[currentKey], matcher)
        .then(result => Object.assign({}, sum, { [currentKey]:result }))
    }, {})
  }


  return matcher(obj)
}

module.exports = applicator
