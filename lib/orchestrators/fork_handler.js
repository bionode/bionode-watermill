'use strict'

const createTask = require('../task.js')

const forkHandler = (task, tasks, breakOut, i, dispatch) => {
  //if (task.info.type === 'fork') {
    console.log('Entered fork!')
    // We now need to duplicate everything after for as many tasks in fork()
    const joinages = []
    task.info.tasks.forEach((forkee, j) => {
      console.log('a forkee: ', forkee.info)
      console.log(forkee.info.type)  //this only when there is other

      const restTasks = tasks.slice(i + 1)
      console.log(restTasks)
      const newRestTasks = restTasks.map(
        t => createTask(dispatch)(
          t.info.props, t.info.operationCreator, hash(t.info.uid + j)
        )
      )

      console.log("fork_bla")
      console.log(forkee)
      console.log(newRestTasks)
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
    reject({type: 'joinage', joinages, currentCtx})   // sketchy
    return breakOut
  //}
}

module.exports = forkHandler