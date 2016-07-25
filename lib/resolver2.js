'use strict'

const Promise = require('bluebird')
const _ = require('lodash')

const resolver = (obj, matcher) => {
  if (_.isArray(obj)) {
    return Promise.all(obj.map(item => resolver(item, matcher)))
  }

  if (_.isString(obj)) {
    return matcher(obj)
  }

  if (_.isPlainObject(obj)) {
    return Promise.reduce(Object.keys(obj), (sum, currentKey) => {
      return resolver(obj[currentKey], matcher)
        .then(result => Object.assign({}, sum, { [currentKey]:result }))
    }, {})
  }
}

module.exports = resolver
