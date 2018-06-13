'use strict'

const PERMANENT_FAILURE = 'PERMANENT_FAILURE'
const TEMPORARY_FAILURE = 'TEMPORARY_FAILURE'
const SUCCESS = 'SUCCESS'

const defaultState = 'starting' // find what conforms to CWL spec..

// Reducer
const reducer = (state = defaultState, action) => {
  switch (action.type) {
    case PERMANENT_FAILURE:
      return 'permanentFailure'
    case TEMPORARY_FAILURE:
      return 'temporaryFailure'
    case SUCCESS:
      return 'success'	 
   default:
      return state
  }
}

// Action creators
reducer.setWorkflowStatusPermanentFailure = () => ({
  type: PERMANENT_FAILURE
})

reducer.setWorkflowStatusTemporaryFailure = () => ({
	type: TEMPORARY_FAILURE
})

reducer.setWorkflowStatusSuccess = () => ({
	type: SUCCESS
})
// Also make action types importable elsewhere
reducer.PERMANENT_FAILURE = PERMANENT_FAILURE

reducer.TEMPORARY_FAILURE = TEMPORARY_FAILURE

reducer.SUCCESS = SUCCESS

module.exports = reducer