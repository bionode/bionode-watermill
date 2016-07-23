'use strict'

// Reducers - this implies there is reducer proxying in this reducer
const taskReducer = require('./task.js')

// Actions
const CREATE_TASK = 'tasks/create'

// Internal Methods

/**
 * updateTask
 * Return a new state with a new task.
 * @param state {Object} reference to current state
 * @param task {Object} the task to add/update into state. Must have `task.uid`.
 */
const updateTask = (state, task) => {
  const { uid } = task
  if (!uid) {
    // TODO test case for this. It is even ever possible for this to happen,
    // if caught properly above?
    throw new Error('task does not have a uid property')
    return state
  }

  return Object.assign({}, state, { [uid]:task })
}

// The `tasks` store item is an object keyed by `task` uids. Each uid points to
// a `task` instance, which also contains `task.uid` for  redundancy. Thus, the
// default state is no tasks.
const defaultState = {}

const reducer = (state = defaultState, action) => {
  const { type } = action

  if (type.startsWith('task/')) {
    // An action has been dispatched for an *assumed existing* task.
    // Proxy the task reducer behind this one.

    // Get the action uid so we can retrieve the current state of that task.
    // If it is not provided, throw and return current state.
    // TODO test case for this
    const { uid } = action
    if (!uid) {
      console.log(action)
      throw new Error('task action not provided with uid')
      return state
    }

    // Use uid to retrieve current task. Throw and return if it is not found.
    // TODO test case for this
    const currentTask = state[uid]
    if (!currentTask) {
      throw new Error('task uid does not match an existing task')
      return state
    }

    // Call taskReducer with current state of corresponding task. This is the
    // point of proxying.
    const newTask = taskReducer(currentTask, action)

    return updateTask(state, newTask)
  }

  switch (type) {
    case CREATE_TASK:
      // We call the task reducer with an undefined state and an empty action.
      // We cannot dispatch this since the task reducer does not have an
      // associated item in the root of the store. This creates a new task with
      // default settings for this waterwheel instance. Once we have this task,
      // which has the `task.uid` property, we can add it to the object of
      // tasks. This is also a proxy to the task reducer as above, yet we need
      // to handle it less directly as the uid is required.
      const { task, cb } = action
      if (!task) {
        throw new Error('task not provided with create task action')
        return state
      }
      if (!cb) {
        throw new Error('callback not provided with create task action')
        return state
      }

      // TODO reconsider. pass action through? what should state be?
      // prefix action with 'proxy'?
      const newTask = taskReducer(task, {})

      // Pass the new task uid back to action creator so that it can resolve
      // with it, and other tasks can use it. The uid is required for every
      // task action.
      cb(newTask.uid)

      return updateTask(state, newTask)
    default:
      return state
  }
}

reducer.createTask = (task) => (dispatch) => new Promise((resolve, reject) => {
  // TODO try/catch this?
  dispatch({
    type: CREATE_TASK,
    task,
    cb: (uid) => resolve(uid)
  })
})

module.exports = reducer
