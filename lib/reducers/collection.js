'use strict'

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
  const { trajectory } = action

  return state
}

reducer.addOutput = (output, trajectory) => (dispatch, getState) => new Promise((resolve, reject) => {
  dispatch({
    type: ADD_OUTPUT,
    trajectory
  })

  resolve()
})

module.exports = reducer
