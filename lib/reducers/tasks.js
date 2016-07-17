'use strict'


const ADD_TASK = 'ADD_TASK'

const defaultState = {} 
const reducer = (state = defaultState, action) => {
  switch (action.type) {
    case ADD_TASK:
      const newTask = {}
      newTask[action.task.uid] = action.task
      return Object.assign({}, state, newTask)
      // return state.concat(taskReducer.createTask(action.task))
      break
    default:
      return state
  }
}

reducer.addTask = (task) => (dispatch) => {
  // Before we can add a task, need to create one so we have its uid
  const { createTask } = require('./task.js')
  dispatch(createTask(task)).then((task) => {
    dispatch({
      type: ADD_TASK,
      task
    }) 
  })
  
  // dispatch({
  //   type: ADD_TASK,
  //   task
  // })
}

module.exports = reducer
