'use strict'

const path = require('path')
const globby = require('globby')
const { bindActionCreators } = require('redux')

const REDUCERS_PATH = './reducers/'

const getActions = (dispatch) => {
  const patterns = ['**/*.js', '!index.js'].map(str => str[0] === '!' ?
    '!' + path.resolve(__dirname, REDUCERS_PATH + str.substring(1)) :
    path.resolve(__dirname, REDUCERS_PATH + str)
  )

  const reducers = globby.sync(patterns)
  
  const actionCreators = reducers.map(reducer => bindActionCreators(Object.assign({}, require(reducer)), dispatch))
  
  const actions = actionCreators.reduce((allActions, actions) => {
    Object.assign(allActions, actions)

    return allActions
  }, {})

  return actions
}

module.exports = getActions
