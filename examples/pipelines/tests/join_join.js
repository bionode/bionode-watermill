'use strict'

// === WATERMILL ===
const {
  task,
  join,
  junction,
  fork
} = require('../../..')

const task0 = task({name: 'task0'}, () => `echo "something0"`)

const task2 = task({name: 'task2'}, () => `echo "something2"`)

const task3 = task({name: 'task3'}, () => `echo "something3"`)

const task4 = task({name: 'task4'}, () => `echo "something4"`)

const task5 = task({name: 'task5'}, () => `echo "something5"`)

const task6 = task({name: 'task6'}, () => `echo "something6"`)

const pipeline = join(
  task0,
  join(
    task5,
    join(
      task4,
      task5 //repetitive task... for test purposes
    ),
    task6
  ),
  task2,
  task3
)

pipeline()