'use strict'

const _ = require('lodash')
const Graph = require('graph.js/dist/graph.full.js')

// Setup server

const http = require('http')
const fs = require('fs')
const path = require('path')

// define path to index.html, otherwise it will be relative to dir where
// scripts is run
const scr = path.join(__dirname, '../../viz/index.html')

// Loading the index file . html displayed to the client

const server = http.createServer(function(req, res) {
  fs.readFile(scr, 'utf-8', function(err, content) {
    res.writeHead(200, {"Content-Type": "text/html"})
    res.end(content)
  })
})

// Loading socket.io

const io = require('socket.io').listen(server)

// sets MaxListeners without limit, important if several instances of
// socket.io are expected (e.g. if we have several junctions, forks, and
// even tasks)
io.sockets.setMaxListeners(0)

// When a client connects, we note it in the console
io.sockets.on('connection', function (socket) {
  console.log('A client is connected!')
})
server.listen(8084)

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
  const graphsonObj = {
    graph: {
      mode: "NORMAL",
      vertices: [],
      edges: []
    }
  }

  // iterate through vertices
  // in ECMAScript 6, you can use a for..of loop
  for (let [key, value] of graph.vertices()) {
    // iterates over all vertices of the graph
    graphsonObj.graph.vertices.push({
      _id: key,
      _type: "vertex",
      values: value
    })
  }
  // iterate through edges
  // in ECMAScript 6, you can use a for..of loop
  for (let [from, to, value] of graph.edges()) {
    // iterates over all vertices of the graph
    graphsonObj.graph.edges.push({
      source: from,
      target: to,
      _type: "edge"
    })
  }

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

  // writes each graphson to file
  fs.writeFile('graphson.json', JSON.stringify(graphsonObj, null, 2), (err) => {
    if (err) throw err
    console.log('File was saved!')
  })

  // Output to visualization server
  io.sockets.on('connection', function (socket) {
    socket.emit('message', graphsonObj)
  })

  Object.keys(obj).forEach(key => Object.assign(obj[key], travelNode(key)))

  return obj
}

function addOutputHandler (state, action) {
  // TODO remove duplicate code b/w this, and creating the next trajection in join
  const {
    type,
    uid,
    resolvedOutput,
    resolvedInput,
    params,
    context,
    name,
    operationString
  } = action
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
  // only adds to graph if uid has in fact a name (from tasks);
  // if resolvedInput doesn't exist and it is outputted to the graph this
  // breaks the pipeline
  if (resolvedInput) {
    // node that in ES6 { type } === { type: type }
    const tempVertex = {
      type,
      name,
      resolvedInput,
      resolvedOutput,
      params,
      operationString //this will pass operationString even if it is undefined
    }
    graph.addNewVertex(uid, tempVertex) //addNewVertex(key, [value])
  } else {
    const tempVertex = { type, name, resolvedOutput, params }
    graph.addNewVertex(uid, tempVertex)
  }

  // Add edge from last trajectory node to current output if no params,
  // otherwise to params (which points to current output)
  if (trajectory.length > 0) {
      graph.addNewEdge(trajectory[trajectory.length - 1], uid)
  }

  return graph
}

function addJunctionVertexHandler(state, action) {
  const { uid, results } = action

  // Clone graph to maintain immutable state
  const graph = state.clone()

  const values = []
  const lastUids = []
  for (const result of results) {
    if (result.tasks) {
      // Dealing with a join/fork/junction
      // Since param nodes never have values, can use these
      result.context.trajectory.forEach((taskUid, i) => {
        values.push(graph.vertexValue(taskUid))
        if (i === result.context.trajectory.length - 1) lastUids.push(taskUid)
      })
    } else {
      // Dealing with a single task
      const taskUid = result.uid
      values.push(graph.vertexValue(taskUid))
      lastUids.push(taskUid)
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
  name: taskState.name,
  resolvedOutput: taskState.resolvedOutput,
  resolvedInput: taskState.resolvedInput,
  context: taskState.context,
  params: taskState.params,
  operationString: taskState.operationString
})

reducer.ADD_JUNCTION_VERTEX = ADD_JUNCTION_VERTEX
reducer.addJunctionVertex = (uid, results) => ({
  type: ADD_JUNCTION_VERTEX,
  uid,
  results
})

reducer.jsonifyGraph = jsonifyGraph

module.exports = reducer
