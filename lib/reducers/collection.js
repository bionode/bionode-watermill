'use strict'

const _ = require('lodash')
const applicator = require('../utils/applicator.js')
const Graph = require('graph.js/dist/graph.full.js')

// Actions
const ADD_OUTPUT = 'collection/add-output'

const defaultState = new Graph()

const reducer = (state = defaultState, action) => {
  switch (action.type) {
    case ADD_OUTPUT:
      return addOutputHandler(state, action)
      break
    default:
      return state
  }
}

function addOutputHandler (state, action) {
  const { uid, trajectory, output } = action

  console.log('trajectory: ', trajectory)
  console.log('output: ', output)

  if (!_.isArray(trajectory)) {
    console.log('Trajectory is not an array!')
    return state
  }

  // TODO add edge from last trajectory node to current output
  const graph = state.clone()
  graph.addNewVertex(uid, output)

  return graph
}

reducer.addOutput = (uid) => (dispatch, getState) => new Promise((resolve, reject) => {
  const { resolvedOutput, trajectory } = getState().tasks[uid]

  dispatch({
    type: ADD_OUTPUT,
    uid,
    output: resolvedOutput,
    trajectory
  })

  resolve()
})

module.exports = reducer
