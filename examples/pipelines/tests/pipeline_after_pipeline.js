'use strict'

// === WATERMILL ===
const {
  task,
  join,
  junction,
  fork
} = require('../../..')

const task0 = task({name: 'task0'}, () => `echo "something0"`)

const task1 = task(
  {
    name: 'task1',
    output: '*.txt'
  }, () => `touch task1_middle.txt`)

const task2 = task(
  {
    name: 'task2',
    output: '*.txt'
  }, () => `touch task2_middle.txt`)

const task3 = task(
  {
    name: 'task3',
    output: '*.txt',
    input: '*.txt'
  }, ({ input }) => `touch ${input}_final`)

const task4 = task(
  {
    name: 'task4',
    input: '*_middle.txt',
    output: '*.txt'
  }, ({ input }) => `echo "something4" >> ${input}`)


const pipeline1 = join(task0, fork(task1, task2), task3)

pipeline1().then(task4)