'use strict'

const _ = require('lodash')
const Promise = require('bluebird')

const { defaultContext } = require('../constants/default-task-state.js')
const { mergeCtx } = require('../ctx')
const hash = require('../utils/hash.js')
const createTask = require('../task.js')

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
        if (task.info.type === 'fork') {
          console.log('Entered fork!')
          // We now need to duplicate everything after for as many tasks in fork()
          const joinages = []
          task.info.tasks.forEach((forkee, j) => {
            console.log('a forkee: ', forkee.info)

            // gets upstream tasks for outer fork
            // if nothing exists after fork this will be an empty array
            const restTasks = tasks.slice(i + 1) // i + 1 retrieves every an
            // array of every task from i to the end of the array, so all
            // tasks after fork will be put into this array (restTasks)

            const newRestTasks = taskCreationDispatcher(dispatch, restTasks, j)

            if (forkee.info.type === 'fork') {
              forkee.info.tasks.forEach((subTask, k) => {
                console.log('fork tasks: ', subTask.info)

                // gets upstream tasks for inner fork and puts them into an
                // array.
                // if nothing exists after fork this will be an empty array
                const upStreamTasks = forkee.info.tasks.slice(k + 1)
                console.log('upStreamTasks: ', upStreamTasks)
                // gets downstream tasks for inner fork and puts them into an
                // array.
                const downStreamTasks = forkee.info.tasks.slice(0, k)
                console.log('downStreamTasks: ', downStreamTasks)

                // executed when there is a fork inside a fork
                if (subTask.info.type === 'fork') {
                  console.log('fork inside fork!')
                  subTask.info.tasks.forEach((subForkee, j2) => {
                    console.log('a subForkee: ', subForkee.info)

                    // gets new uid for tasks outside first fork whenever
                    // they encounter a inner fork and so they must have to
                    // levels of uid modification (j and j2, for outer and
                    // inner fork respectively)
                    const outerRestTasks = taskCreationDispatcher(dispatch, newRestTasks, j + j2)

                    // creates new uid for upstream tasks of the inner fork
                    // e.g. (task0, fork(task1, task2) -> upstream task is
                    // an array with task0
                    const newUpTasks = taskCreationDispatcher(dispatch, upStreamTasks, j2)

                    // creates a new uid for downstream tasks of the inner fork
                    // e.g. (task0, fork(task1, task2), task3) -> downstream
                    // task is an array with task3
                    const newDownTasks = taskCreationDispatcher(dispatch, downStreamTasks, j2)

                    // concat downstream tasks of inner fork inside fork +
                    // inner fork + upstream tasks of inner fork +
                    // upstream tasks of outer fork (newRestTasks)
                    const lineage = newDownTasks.concat([subForkee]).concat(newUpTasks).concat(outerRestTasks)
                    lineageGetter(dispatch, lineage, joinages)
                  })
                } else if (subTask.info.type !== undefined) {
                  // executed when there is a junction or join inside a fork
                  console.log(subTask.info.type + ' inside fork!')
                  const lineage = downStreamTasks.concat([subTask]).concat(upStreamTasks).concat(newRestTasks)
                  lineageGetter(dispatch, lineage, joinages)
                }
              })
            } else {
              console.log('simple task!')
              // performed whenever there is no handling of arrays of tasks
              // within fork orchestrator (join, junction ,fork)

              // concats fork task with tasks after fork, when fork is
              // simple (fork(task1, task2)
              const lineage = [forkee].concat(newRestTasks)
              lineageGetter(dispatch, lineage, joinages)
            }

          })
          breakOut = true
          reject({ type: 'joinage', joinages, currentCtx })   // sketchy
          return
        }
      }

      // Call next task with currentCtx
      if (!breakOut) {
        task(_.noop, currentCtx).then((results) => {
          console.log("results: ", results)
          // Resolve to a new context with a ctx merge strategy
          let newCtx
          if (!_.isArray(results)) {
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
      tasks: tasks.map(task => task),  //returns an array of tasks
      type: 'join'  // used to check the type of orchestrator
    }

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
const taskCreationDispatcher = (dispatch, taskInstance, number) => {
  const taskNewUid = taskInstance.map(
    t => createTask(dispatch)(
      t.info.props, t.info.operationCreator, hash(t.info.uid + number)
    )
  )
  return taskNewUid
}


module.exports = join
