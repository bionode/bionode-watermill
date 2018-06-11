'use strict'

const PERMANENT_FAILURE = 'PERMANENT_FAILURE'

const defaultState = 'starting' // find what conforms to CWL spec..

// Reducer
const reducer = (state = defaultState, action) => {
  switch (action.type) {
    case PERMANENT_FAILURE:
      return 'permanentFailure'
    case SOME_ACTION:
  	  return someActionHandler(state, action)
  	  break
   default:
      return state
  }
}

// Action creators
reducer.setWorkflowStatusPermanentFailure = () => ({
  type: PERMANENT_FAILURE
})
// Also make action types importable elsewhere
reducer.PERMANENT_FAILURE = 'PERMANENT_FAILURE'

module.exports = reducer