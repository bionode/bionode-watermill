# Orchestration

A **task** is the fundamental **unit** for building pipelines. It

* has a **unique input/output pair** defined by glob pattern(s) and/or streams
* has a **single params object**. Params should be used when they can be applied
  to a task with the same I/O but alter the output. e.g. command flags
* **emits** `task.finish` with a reference to that task object in state

*Example*

```javascript
/*
 * Usage: samtools index [-bc] [-m INT] <in.bam> [out.index]
 * Options:
 *   -b       Generate BAI-format index for BAM files [default]
 *   -c       Generate CSI-format index for BAM files
 *   -m INT   Set minimum interval size for CSI indices to 2^INT [14]
 */

const samtoolsIndex = task({
  input: '*.bam',
  output: '*.bai',
  params: { format: 'b' }
  name: 'samtools index'
}, ({ input, params }) => shell(`samtools index -${params.format} ${input}`))

samtoolsIndex()
  .on('task.finish', (results) => console.log(results.resolvedOutput))
```

A **join** is used to run a **sequence of tasks in order**. It

* must pass the `trajectory` between tasks (including the trajectory passed into `join`)
* **emits** `join.finish` with a `trajectory`
* returns/emits reference to each task

*Example*

```javascript
// task1 will resolve input to filesystem,
// task2 will resolve input to files at collection nodes defined by its trajectory
// (which will be created based on the trajectory task1 returns and join passes to task2)
const joined = join(task1, task2)

joined()
  .tasks((t1, t2) => {
    t1.on('task.finish', (results) => {
      anotherTask(results.trajectory)()
        .on('task.finish', () => console.log('another task finished'))
    })
  })
  .on('join.finish', (results) => console.log(results.trajectory))
// Which looks like:
// task1 -> task2
// task1 -> anotherTask

// NOTE: This is only recommended for advanced or dynamic cases. Try to avoid
// retreiving task references from join or parallel.
// Another way to do it:
const pipeline = join(task1, parallel(task2, anotherTask), task3)
// However, in this case, task3 will only start after task2 and anotherTask finish
// task1 -> task2        ⌉
//                       | -> task3
// task1 -> anotherTask  ⌋
// Alternatively you could do:
const pipeline = join(task1, parallel(join(task2, task3), anotherTask))
// task1 -> task2 -> task3
// task1 -> anotherTask
```

A **parallel** is used to run **a set of tasks simultaneously**. It

* must pass the trajectory passed into `parallel` into **each task**
* **emits** `parallel.finish` with a list of trajections
* returns/emits reference to each task

```javascript
const parallelized = parallel(taskA, taskB)

parallelized()
  .tasks((a, b) => { /* do stuff if u want */ })
  .on('parallel.finish', (results) => {
    console.log(results.trajectory) // an array
  })
```

