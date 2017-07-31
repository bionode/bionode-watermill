'use strict'

// === WATERMILL ===
const {
  task,
  join,
  junction,
  fork
} = require('../../..')

const task0 = task({name: 'coco'}, () => `echo "something0"`)

const task1 = task({name: 'xixi'}, () => `echo "something1"`)

const task2 = task({name: 'foo'}, () => `echo "something2"`)

const task3 = task({name: 'bar'}, () => `echo "something3"`)

const task4 = task({name: 'test'}, () => `echo "something4"`)

const task5 = task({name: 'test1'}, () => `echo "something5"`)

const task6 = task({name: 'test6'}, () => `echo "something6"`)

const pipeline = join(task0, join(task1, join(task4, task5), task6), task2, task3)

pipeline()