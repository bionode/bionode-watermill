'use strict'

// === WATERMILL ===
const {
  task,
  join,
  junction,
  fork
} = require('../../..')

const task0 = task({name: 'coco'}, () => `echo "something0"`)

const task1 = task({name: 'nut'}, () => `echo "something1"`)

const task2 = task({name: 'foo'}, () => `echo "something2"`)

const task3 = task({name: 'bar'}, () => `echo "something3"`)

const task4 = task({name: 'test'}, () => `echo "something4"`)

const task5 = task({name: 'test1'}, () => `echo "something5"`)

const pipeline = join(task0, fork(join(task4, fork(task1, task3)), task2), task5)

pipeline()