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

const jsonifyGraph = (graph, obj = {}) => {
  const assignVal = (obj, key, val) => {
    if (val === undefined) {
      obj[key] = { _output: '' }
    } else {
      obj[key] = { _output: val }
    }
  }

  for (const [key, val] of graph.sources()) {
    assignVal(obj, key, val)
  }

  const travelNode = (node) => {
    const obj = {}
    for (const [key, val] of graph.verticesFrom(node)) {
      assignVal(obj, key, val)
      Object.assign(obj[key], travelNode(key))
    }

    return obj
  }

  Object.keys(obj).forEach(key => {
    Object.assign(obj[key], travelNode(key))
  })

  return obj
}

function addOutputHandler (state, action) {
  const { uid, trajectory, output } = action

  // console.log(`trajectory for ${uid}: `, trajectory)
  // console.log('output: ', output)

  if (!_.isArray(trajectory)) {
    console.log('Trajectory is not an array!')
    return state
  }

  // TODO add edge from last trajectory node to current output
  const graph = state.clone()
  graph.addNewVertex(uid, output)

  if (trajectory.length > 0) {
    graph.addNewEdge(trajectory[trajectory.length - 1], uid)
  }

  console.log('collection as json: ', JSON.stringify(jsonifyGraph(graph), null, 2))

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

  resolve(uid)
})

module.exports = reducer
