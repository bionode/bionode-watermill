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

      // TODO better checking if task/join/junction/fork
      if (task.info) {
        //if (task.info.props) {
        //  console.log('Joining to task: ' + task.info.props.name)
        //}
        console.log('Joining to task: ' + _.get(task, 'info.props.name'))
        // Add this task to list of uids for this join
        uids.push(task.info.uid)
        if (task.info.type === 'fork') {
          // console.log('ABOUT TO FORK, currentCtx is: ', currentCtx)
          // We now need to duplicate everything after for as many tasks in fork()
          const joinages = []
          task.info.tasks.forEach((forkee, j) => {
            console.log('a forkee: ', forkee.info)
            const restTasks = tasks.slice(i+1)
            // console.log(restTasks[0].info.uid)
            const newRestTasks = restTasks.map(t => createTask(dispatch)(t.info.props, t.info.operationCreator, hash(t.info.uid + j)))
            // restTasks.forEach(t => t.setInstanceUid(hash(t.info.uid + j)))
            // console.log('restTasks: ', newRestTasks.map(t => t.info.uid))
            // Join down the line
            // Ideally inside here, other forks will also be handled
            const lineage = [forkee].concat(newRestTasks)
            console.log('lineage: ', lineage.map(t => t.info.uid))
            const joinage = join(dispatch).apply(null, lineage)
            joinages.push(joinage)
          })
          breakOut = true
          // Run each new joinage
          // joinages.forEach(joinage => joinage(currentCtx))
          // Promise.all(joinages.map(j => j(currentCtx))).then((results) => {
          //   console.log('results: ', results)
          // })
          reject({ type: 'joinage', joinages, currentCtx })
          return
        }
      }

      // Call next task with currentCtx
      if (!breakOut) {
        task(_.noop, currentCtx).then((results) => {
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

    joinInvocator.info = tasks.map(task => task.info)

    return joinInvocator
  }
}

module.exports = join
