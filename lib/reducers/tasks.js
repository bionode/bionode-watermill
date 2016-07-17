'use strict'

const ADD_TASK = 'ADD_TASK'

const reducer = (state = [], action) => {
  switch (action.type) {
    case ADD_TASK:
      const { task } = action

      return state.concat(task)
      break
    default:
      return state
  }
}

reducer.addTask = (task) => ({
  type: ADD_TASK,
  task
})

module.exports = reducer
