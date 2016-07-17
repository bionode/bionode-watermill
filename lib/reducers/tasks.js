'use strict'

const { _createTask } = require('./task.js')

const ADD_TASK = 'ADD_TASK'

const defaultState = {} 
const reducer = (state = defaultState, action) => {
  switch (action.type) {
    case ADD_TASK:
      // Secretly "dispatches" create/task to task reducer
      const newTask = _createTask(action.task)
      const obj = {}
      obj[newTask.uid] = newTask
      return Object.assign({}, state, obj)
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
