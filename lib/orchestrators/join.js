'use strict'

const _ = require('lodash')
const Promise = require('bluebird')

const { defaultContext } = require('../constants/default-task-state.js')
const { mergeCtx } = require('../ctx')
const hash = require('../utils/hash.js')
const createTask = require('../task.js')
const orchestratorUidGenerator = require('./orchestrator-uid-generator')

// a variable that should be global in order to store outermostTasks of fork
// which should be wrapped in a join for now, since we need a way to fetch
// the tasks outside the scope of fork and join is useful to control the
// flow before fork
let outermostTasks

// this variable allows to check whether it is the firstFork inside join or
// if there is some level of fork inside other forks.
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
        //console.log('taskuidsjoin: ', uids)
      }

      if (task.info.type === 'fork') {
        //console.log('fork encountered inside join! ', task.info, firstFork)
        const joinages = []
        // We now need to duplicate everything after for as many tasks in fork()
        task.info.tasks.forEach((forkee, j) => {
          //console.log('a forkee: ', forkee.info)
          // gets downStream tasks for inner fork and puts them into an
          // array. if nothing exists after fork this will be an empty array
          const downStreamTasks = tasks.slice(i + 1)
          // makes an hash for upStreamTasks
          const upStreamTasks = tasks.slice(0, i)
          // puts the array of uids into a unique uid for each downStremTasks
          const upStreamUids = upStreamTasks.map(
            t => t.info.uid
          ).reduce((a, b) => hash(a + b), 0)
          // array of every task from i to the end of the array, so all
          // tasks after fork will be put into this array (restTasks)
          if (firstFork === true) {
            // this if statement is used to store outermostTasks for inner
            // orchestrators inside fork
            const newdownStreamTasks = taskCreationDispatcher(dispatch, downStreamTasks, task.info.uid, forkee.info.uid, upStreamUids)
            // if a task is found for outer fork it will dispatch the outer
            // task as newdownStreamTasks, otherwise it will concat the other
            // orchestrators
            const lineage = [forkee].concat(newdownStreamTasks)
            lineageGetter(dispatch, lineage, joinages)
            // this is used to control when firstFork instance ends
            if (j === task.info.tasks.length - 1) {
              firstFork = false
              outermostTasks = newdownStreamTasks
            }
          } else {
            if (forkee.info.props) {
              // this handles if fork is within fork without join wrapping
              // the inner fork
              const sumOuterMostTasks = outermostTasks.map(
                t => t.info.operationCreator
              ).reduce((a, b) => hash(a + b), 0)

              // creates an hash for all downStreamTasks
              const sumDownStreamTasks = downStreamTasks.map(
                t => t.info.operationCreator    // it uses task operationCreator
              ).reduce((a, b) => hash(a + b), 0)
              // creates an hash for all outermostTasks
              //console.log('taskstotest:',
              // sumOuterMostTasks,sumDownStreamTasks)
              // then checks if the hash is the same
              // it will be the same if fork is followed by another fork
              // it will be different if fork is wrapped by a join inside a fork
              if (sumOuterMostTasks !== sumDownStreamTasks) {
                const newdownStreamTasks = taskCreationDispatcher(dispatch, downStreamTasks, task.info.uid, forkee.info.uid, upStreamUids)
                const newOutermostTasks = taskCreationDispatcher(dispatch, outermostTasks, task.info.uid, forkee.info.uid, upStreamUids)
                // notice that if not within a fork it creates a lineage
                // with the joint tasks plus the outermostTask with a new uid
                const lineage = [forkee].concat(newdownStreamTasks).concat(newOutermostTasks)
                lineageGetter(dispatch, lineage, joinages)
              } else { // if a fork is found within another fork
                // this avoids the duplication of the last task
                const newOutermostTasks = taskCreationDispatcher(dispatch, outermostTasks, task.info.uid, forkee.info.uid, upStreamUids)
                // in the case of fork it simply adds to the forkee task the
                // outermost task
                const lineage = [forkee].concat(newOutermostTasks)
                lineageGetter(dispatch, lineage, joinages)
              }
            } else if (forkee.info.type === 'fork') {
              // handles everything that is not task out of the scope of
              // firstFork instance
              // this handles if fork is within fork without join wrapping
              // the inner fork
              const sumOuterMostTasks = outermostTasks.map(
                t => t.info.operationCreator
              ).reduce((a, b) => hash(a + b), 0)

              // creates an hash for all downStreamTasks
              const sumDownStreamTasks = downStreamTasks.map(
                t => t.info.operationCreator    // it uses task operationCreator
              ).reduce((a, b) => hash(a + b), 0)
              // creates an hash for all outermostTasks
              //console.log('taskstotest:',
              // sumOuterMostTasks,sumDownStreamTasks)
              // then checks if the hash is the same
              // it will be the same if fork is followed by another fork

              if (sumOuterMostTasks !== sumDownStreamTasks) {
                const newdownStreamTasks = taskCreationDispatcher(dispatch, downStreamTasks, task.info.uid, forkee.info.uid, upStreamUids)
                const newOutermostTasks = taskCreationDispatcher(dispatch, outermostTasks, task.info.uid, forkee.info.uid, upStreamUids)
                // notice that if not within a fork it creates a lineage
                // with the joint tasks plus the outermostTask with a new uid
                const lineage = [forkee].concat(newdownStreamTasks).concat(newOutermostTasks)
                lineageGetter(dispatch, lineage, joinages)
              } else { // if a fork is found within another fork
                // this avoids the duplication of the last task
                const newOutermostTasks = taskCreationDispatcher(dispatch, outermostTasks, task.info.uid, forkee.info.uid, upStreamUids)
                // in the case of fork it simply adds to the forkee task the
                // outermost task
                const lineage = [forkee].concat(newOutermostTasks)
                lineageGetter(dispatch, lineage, joinages)
              }
            } else {
              // for join and junction immediately after a fork in order to
              // not duplicate last task
              const newOutermostTasks = taskCreationDispatcher(dispatch, outermostTasks, task.info.uid, forkee.info.uid, upStreamUids)
              const lineage = [forkee].concat(newOutermostTasks)
              lineageGetter(dispatch, lineage, joinages)
            }
          }
        })
        breakOut = true
        reject({ type: 'joinage', joinages, currentCtx })
        return
      }

      // Call next task with currentCtx
      if (!breakOut) {
        task(_.noop, currentCtx).then((results) => {
          //console.log("results: ", results)
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
  //console.log('lineage: ', lineage.map(t => t.info.uid))
  const joinage = join(dispatch).apply(null, lineage)
  joinages.push(joinage)
  //console.log('joinages: ', joinages)
}

// this function generates a new uid and helps duplicate task for fork
// this basically takes a taskInstance and a number coming from a loop in
// the fork's task
const taskCreationDispatcher = (dispatch, taskInstance, forkUid, forkeeUid, upStreamUid) => {
  //console.log("taskInstance: ", taskInstance, forkUid, forkeeUid, upStreamUid)
  const taskNewUid = taskInstance.map(
    t => createTask(dispatch)(
       t.info.props, t.info.operationCreator, hash(t.info.uid + forkUid + forkeeUid + upStreamUid)
    )
  )
  //console.log('taskNewUid: ', taskNewUid)
  return taskNewUid
}


module.exports = join
