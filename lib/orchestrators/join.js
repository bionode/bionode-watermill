'use strict'

const _ = require('lodash')
const Promise = require('bluebird')

const { defaultContext } = require('../constants/default-task-state.js')
const { mergeCtx } = require('../ctx')
const hash = require('../utils/hash.js')
const createTask = require('../task.js')

const join = (dispatch, logger) => {
  return function (...tasks) {
    console.log("join task")
    console.log(tasks)
    let uids = []
    let breakOut = false
    const accumulator = (currentCtx, task, i, length) => new Promise((resolve, reject) => {
      if (breakOut) return

      // TODO better checking if task/join/junction/fork
      if (task.info) {
        console.log("check")
        // console.log(task) // this retrieves the function
        console.log(task.info) // this retrieves an object for tasks and an
        // array for orchestrators
        //console.log(task.info.type)

        // enables checking of orchestrator versus task
        if (_.isArray(task.info.tasks)) {
          console.log(task.info.type +' found! Tasks array: ', task.info.tasks)
          console.log(task)
        } else {
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
            const restTasks = tasks.slice(i + 1)
            console.log(restTasks)


            console.log("fork_bla")
            console.log(tasks)
            console.log(forkee)

            if (forkee.info.tasks) {
              forkee.info.tasks.forEach((subTask, k) => {
                console.log('fork tasks: ', subTask.info)
                console.log(k)
                console.log("fork_bla2")
                console.log(forkee.info.tasks)

                if (subTask.info.type === 'fork') {
                  console.log("yes_array")
                  console.log(subTask.info.tasks)

                  // gets upstream tasks for inner fork
                  // if nothing exists after fork this will be an empty array
                  const restSubTasks = forkee.info.tasks.slice(k + 1)
                  console.log(restSubTasks)


                  subTask.info.tasks.forEach((subForkee, j2) => {
                    console.log('a subForkee: ', subForkee.info)
                    // we can use

                    const newRestTasks = restTasks.map(
                      t => createTask(dispatch)(
                        t.info.props, t.info.operationCreator, hash(t.info.uid + j + j2)
                      )
                    )

                    console.log(subForkee)
                    const lineage = [subForkee].concat(newRestTasks)
                    console.log('lineage: ', lineage.map(t => t.info.uid + j2))
                    const joinage = join(dispatch).apply(null, lineage)
                    joinages.push(joinage)
                    console.log('joinages :', joinages)
                  })
                }
              })
            } else {
              const newRestTasks = restTasks.map(
                t => createTask(dispatch)(
                  t.info.props, t.info.operationCreator, hash(t.info.uid + j)
                )
              )
              console.log("no_array")
              console.log(forkee)
              const lineage = [forkee].concat(newRestTasks)
              console.log('lineage: ', lineage.map(t => t.info.uid))
              const joinage = join(dispatch).apply(null, lineage)
              joinages.push(joinage)
              console.log('joinages :', joinages)
            }

            // const lineage = [forkee].concat(newRestTasks)
            // console.log('lineage: ', lineage.map(t => t.info.uid))
            // const joinage = join(dispatch).apply(null, lineage)
            // joinages.push(joinage)
          })
          breakOut = true

          // Run each new joinage
          // joinages.forEach(joinage => joinage(currentCtx))
          // Promise.all(joinages.map(j => j(currentCtx))).then((results) => {
          //   console.log('results: ', results)
          // })
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

module.exports = join
