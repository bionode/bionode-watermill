'use strict'

const _ = require('lodash')
const applicator = require('../utils/applicator.js')

// Actions
const ADD_OUTPUT = 'collection/add-output'

const defaultState = {}

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

  if (trajectory.length === 0) {
    // Can append to root of state
    const newState = Object.assign({}, state, { [uid]:output })
    console.log('newState: ', newState)
    return newState
  } else {
    // Append to leaf of trajection
  }

  return state
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
