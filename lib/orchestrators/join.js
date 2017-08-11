'use strict'

const _ = require('lodash')
const Promise = require('bluebird')

const { defaultContext } = require('../constants/default-task-state.js')
const { mergeCtx } = require('../ctx')
const hash = require('../utils/hash.js')
const createTask = require('../task.js')
const orchestratorUidGenerator = require('./orchestrator-uid-generator')

let outermostTasks

let firstFork = true

const join = (dispatch, logger) => {
  return function (...tasks) {
    let uids = []
    let breakOut = false

    const accumulator = (currentCtx, task, i, length) => new Promise((resolve, reject) => {
      if (breakOut) return

      if (task.info) {
        // enables checking of orchestrator versus task
        if (_.isArray(task.info.tasks)) {
          // executed for orchestrators
          console.log(task.info.type +' found! Tasks array: ', task.info.tasks)
        } else {
          // executed for tasks
          console.log('Joining to task: ' + _.get(task, 'info.props.name'))
        }
        // Add this task to list of uids for this join
        uids.push(task.info.uid)
        console.log('taskuidsjoin: ', uids)
      }

      if (task.info.type === 'fork') {
        console.log('fork encountered inside join! ', task.info, firstFork)
        const joinages = []
        // We now need to duplicate everything after for as many tasks in fork()
        task.info.tasks.forEach((forkee, j) => {
          console.log('a forkee: ', forkee.info)
          // gets upstream tasks for inner fork and puts them into an
          // array. if nothing exists after fork this will be an empty array
          const upStreamTasks = tasks.slice(i + 1)
          console.log(upStreamTasks)
          // array of every task from i to the end of the array, so all
          // tasks after fork will be put into this array (restTasks)

          if (firstFork === true) {
            console.log('test:', forkee.info, outermostTasks)
            const newUpStreamTasks = taskCreationDispatcher(dispatch, upStreamTasks, task.info.uid, forkee.info.uid)
            const lineage = [forkee].concat(newUpStreamTasks)
            lineageGetter(dispatch, lineage, joinages)
            outermostTasks = newUpStreamTasks
            if (j == task.info.tasks.length - 1) {
              firstFork = false
              outermostTasks = newUpStreamTasks
              console.log(firstFork, newUpStreamTasks)
            }
          } else {
            if (forkee.info.props) {
              console.log('test2:', forkee.info, outermostTasks)
              const newUpStreamTasks = taskCreationDispatcher(dispatch, upStreamTasks, task.info.uid, forkee.info.uid)
              const newOutermostTasks = taskCreationDispatcher(dispatch, outermostTasks, task.info.uid, forkee.info.uid)
              const lineage = [forkee].concat(newUpStreamTasks).concat(newOutermostTasks)
              lineageGetter(dispatch, lineage, joinages)
            }
          }
        })
        breakOut = true
        reject({ type: 'joinage', joinages, currentCtx })   // sketchy
        return
      }

      // Call next task with currentCtx
      if (!breakOut) {
        task(_.noop, currentCtx).then((results) => {
          console.log("results: ", results)
          // Resolve to a new context with a ctx merge strategy
          let newCtx
          if (!_.isArray(results) && !_.isUndefined(results.context.trajectory)) {
            newCtx = mergeCtx('join')(currentCtx, results)
            resolve(newCtx)
          } else {
            resolve({ context: currentCtx, results })
          }
        })
      }
    })

    const joinInvocator = (cb = _.noop, ctx = defaultContext) =>
      Promise.reduce(tasks, accumulator, ctx)
        .then(results => Promise.resolve(Object.assign({}, {
          type: results.type,
          uid: hash(uids.join('')),
          name: 'generic join name',
          tasks: uids,
          context: {
            trajectory: results.trajectory
          }
        })))
        .catch((err) => {
          if (err.type === 'joinage') {
            return Promise.all(err.joinages.map(j => j(() => {}, err.currentCtx)))
          } else {
            throw err
          }
        })
        .asCallback(cb)

    joinInvocator.info = {
      tasks: tasks,  //returns an array of tasks
      type: 'join',  // used to check the type of orchestrator
      uid: hash(orchestratorUidGenerator(tasks).join(''))
    }

    //console.log("joinInvocator")
    //console.log(joinInvocator.info)

    return joinInvocator
  }
}

// function to get a lineage and dispatch it as a join for every fork instance,
// wither they are tasks, join, junction or even other fork
const lineageGetter = (dispatch, lineage, joinages) => {
  console.log('lineage: ', lineage.map(t => t.info.uid))
  const joinage = join(dispatch).apply(null, lineage)
  joinages.push(joinage)
  console.log('joinages: ', joinages)
}

// this function generates a new uid and helps duplicate task for fork
// this basically takes a taskInstance and a number coming from a loop in
// the fork's task
const taskCreationDispatcher = (dispatch, taskInstance, forkUid, forkeeUid) => {
  console.log("taskInstance: ", taskInstance, forkUid, forkeeUid)
  const taskNewUid = taskInstance.map(
    t => createTask(dispatch)(
      t.info.props, t.info.operationCreator, hash(t.info.uid + forkUid + forkeeUid) // somehow outer task is not unique
    )
  )
  return taskNewUid
}


module.exports = join
