'use strict'

// === WATERMILL ===
const {
  task,
  join,
  junction,
  fork
} = require('../../..')

const task0 = task({name: 'task0'}, () => `echo "something0"`)

const task1 = task({name: 'task1'}, () => `echo "something1"`)

const task2 = task({name: 'task2'}, () => `echo "something2"`)

const task3 = task({name: 'task3'}, () => `echo "something3"`)

const task4 = task({name: 'task4'}, () => `echo "something4"`)

const task5 = task({name: 'task5'}, () => `echo "something5"`)

const task6 = task({name: 'task6'}, () => `echo "something6"`)

const task7 = task({name: 'task7'}, () => `echo "something7"`)

const pipeline = join(
  task0,
  fork(
    task4,
    fork(
      task1,
      fork(
        task3,
        fork(task6, task7)
      )
    ),
    task2
  ),
  task5
)

pipeline()