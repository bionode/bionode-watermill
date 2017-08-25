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

const task3_2 = task({name: 'task3_2'}, () => `echo "something3_2"`)

const task8 = task({name: 'task8'}, () => `echo "something8"`)


// alternative with manual handling of task5 (the last leave of each branch)
const pipeline = join(
  task0,
  fork(
    join(task2,task5),
    join(task4,
      fork(
        join(task1,task5),
        join(task3,task5)
      )
    )
  )
)

const pipeline2 = join(
  task0,
  fork(
    join(
      task4,
        fork(
          task1,
          task3
        ),
      task6
    ),
    task2
  ),
  task5
)

const pipeline3 = join(
  task0,
  fork(
      fork(
        task1,
        task4
      ),
      task6
  ),
  task5
)

const pipeline4 = join(
  task0,
  fork(
    join(
      task4,
      fork(
        join(
          task7,
          fork(
            task1,
            task3
          )
        ),
        join(task8,task6)
      )
    ),
    task2
  ),
  task5
)


const pipeline5 = join(
  task0,
  fork(
    join(
      task4,
      fork(
        join(
          task7,
          junction(
            task1,
            task3
          ),
          task6
        ),
        task6
      )
    ),
    task2
  ),
  task5, task3_2
)

const pipeline6 = join(
  task0,
  fork(
    join(
      task4,
      fork(
        join(
          task7,
          fork(
            join(
              task1,
              fork(task3, task3_2)
            ),
            task8
          )
        ),
        task6
      )
    ),
    task2
  ),
  task5
)

// not working
const pipeline7 = join(
  task0,
  fork(task4, task3),
  task5,
  fork(task1, task2),
  task6
)

// working
const pipeline7_2 = join(
  task0,
  fork(join(task4,task5,fork(task1,task2)), join(task3,task5,fork(task1,task2))),
  task6
)

// not working
const pipeline8 = join(
  task0,
  fork(task4, task3),
  task5,
  junction(task1, task2),
  task6
)

// working
const pipeline8_2 = join(
  task0,
  fork(join(task4, task5, junction(task1, task2)), join(task3,task5,junction(task1, task2))),
  //task5,
  //junction(task1, task2),
  task6
)

// edit this line to run the desired pipeline.
// documentation on the usage of these pipelines may be found in the link below
// https://github.com/bionode/GSoC17/blob/master/Journal/Week_11.md
pipeline7_2().then((results) => {
  console.log('RESULTSNEW: ', results)
})