'use strict'

const _ = require('lodash')
const applicator = require('../utils/applicator.js')
const Graph = require('graph.js/dist/graph.full.js')

// Actions
const ADD_OUTPUT = 'collection/add-output'
const ADD_JUNCTION_VERTEX = 'collection/add-junction-vertex'

const defaultState = new Graph()

const reducer = (state = defaultState, action) => {
  switch (action.type) {
    case ADD_OUTPUT:
      return addOutputHandler(state, action)
      break
    case ADD_JUNCTION_VERTEX:
      return addJunctionVertexHandler(state, action)
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
  // TODO remove duplicate code b/w this, and creating the next trajection in join
  const { uid, output, params, context } = action
  const { trajectory } = context

  // console.log(`trajectory for ${uid}: `, trajectory)
  // console.log(`output for ${uid}: `, output)
  // console.log(`params for ${uid}: `, params)

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

  // console.log('collection as json: ', JSON.stringify(jsonifyGraph(graph), null, 2))

  return graph
}

function addJunctionVertexHandler(state, action) {
  const { uid, results } = action

  console.log('gonna do something with: ')
  console.log(JSON.stringify(results, null, 2))

  console.log('gonna make a node with uid: ', uid)

  // Clone graph to maintain immutable state
  const graph = state.clone()

  const values = []
  const lastUids = []
  for (const result of results) {
    if (result.tasks) {
      // Dealing with a join/fork/junction
      // Since param nodes never have values, can use these
      result.tasks.forEach((taskUid, i) => {
        values.push(graph.vertexValue(taskUid))
        if (i === result.tasks.length - 1) lastUids.push(taskUid)
      })
    }
  }
  // TODO flatten values?
  // I think each vertex value is already an array only
  graph.addNewVertex(uid, _.flatten(values))

  // Edge from last task of each to new junction node
  lastUids.forEach(lastUid => graph.addNewEdge(lastUid, uid))

  return graph
}

reducer.ADD_OUTPUT = ADD_OUTPUT
reducer.addOutput = (uid, taskState) => ({
  type: ADD_OUTPUT,
  uid,
  output: taskState.resolvedOutput,
  context: taskState.context,
  params: taskState.params
})

reducer.ADD_JUNCTION_VERTEX = ADD_JUNCTION_VERTEX
reducer.addJunctionVertex = (uid, results) => ({
  type: ADD_JUNCTION_VERTEX,
  uid,
  results
})

reducer.jsonifyGraph = jsonifyGraph

module.exports = reducer
