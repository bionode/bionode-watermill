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
      obj[key] = { _value: '' }
    } else {
      obj[key] = { _value: val }
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

  Object.keys(obj).forEach(key => Object.assign(obj[key], travelNode(key)))

  return obj
}

function addOutputHandler (state, action) {
  const { uid, trajectory, output, params } = action

  console.log(`trajectory for ${uid}: `, trajectory)
  console.log('output: ', output)
  console.log('params: ', params)

  if (!_.isArray(trajectory)) {
    console.log('Trajectory is not an array!')
    return state
  }

  // Clone graph to maintain immutable state
  const graph = state.clone()

  // Add a node for current task
  graph.addNewVertex(uid, output)

  // If params is not {}, add a node for it
  const paramsString = JSON.stringify(params)
  if (paramsString !== '{}') {
    // TODO check if params already exists
    graph.addNewVertex(paramsString)

    // Add an edge from params to current task as well
    graph.addNewEdge(paramsString, uid)
  }

  // Add edge from last trajectory node to current output if no params,
  // otherwise to params (which points to current output)
  if (trajectory.length > 0) {
    if (paramsString === '{}') {
      graph.addNewEdge(trajectory[trajectory.length - 1], uid)
    } else {
      graph.addNewEdge(trajectory[trajectory.length - 1], paramsString)
    }
  }

  console.log('collection as json: ', JSON.stringify(jsonifyGraph(graph), null, 2))

  return graph
}

reducer.addOutput = (uid) => (dispatch, getState) => new Promise((resolve, reject) => {
  const { resolvedOutput, trajectory, params } = getState().tasks[uid]

  dispatch({
    type: ADD_OUTPUT,
    uid,
    output: resolvedOutput,
    trajectory,
    params
  })

  resolve(uid)
})

module.exports = reducer
