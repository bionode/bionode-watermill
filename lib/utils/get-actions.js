'use strict'

const path = require('path')
const globby = require('globby')
const { bindActionCreators } = require('redux')

const REDUCERS_PATH = path.resolve(__dirname, '../reducers') + '/'

const getActions = (dispatch) => {
  const patterns = ['**/*.js', '!index.js'].map(str => str[0] === '!' ?
    '!' + REDUCERS_PATH + str.substring(1) :
    REDUCERS_PATH + str
  )

  // Only happens once, at instantiation of a waterwheel pipeline, so OK to use sync
  const reducers = globby.sync(patterns)

  const actionCreators = reducers.map(reducer => bindActionCreators(Object.assign({}, require(reducer)), dispatch))

  const actions = actionCreators.reduce((allActions, actions) => Object.assign({}, allActions, actions), {})

  return actions
}

module.exports = getActions
