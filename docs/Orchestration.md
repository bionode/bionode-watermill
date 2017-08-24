# Orchestration

Bionode-watermill has three main orchestrators: **join**, **junction** and 
**fork**. In current API, these orchestrators gather a group of tasks that 
are described inside them and apply different effects on the pipeline, in 
which they are used.

## Join

A `join` is used to run a **sequence of tasks in order**. It

* must pass the `trajectory` between tasks (including the trajectory passed into `join`)
* **emits** `join.finish` with a `trajectory`
* returns/emits reference to each task

*Example*

```javascript
// task1 will resolve input to filesystem,
// task2 will resolve input to files at collection nodes defined by its trajectory
// (which will be created based on the trajectory task1 returns and join passes to task2)
const pipeline = join(task1, task2)
// executes the pipeline
pipeline()
```

The above represented `pipeline` will first run `task1` and **only after 
`task1`** is
 finished, starts running `task2`.

## Junction

A ``junction`` is used to run **a set of tasks simultaneously** and **waits 
for the results** of all tasks within junction before running anythin else.
. It

* must pass the trajectory passed into `junction` into **each task**
* **emits** `junction.finish` with a list of trajections
* returns/emits reference to each task

*Example*

```javascript
// executes both taskA and taskB at the same time
// however one task may finish before the other depending on the tasks itself
const pipeline2 = junction(taskA, taskB)
// executes the pipeline
pipeline2()
```

The above represented `pipeline2` will run both tasks (`taskA` and `taskB`) at the 
same time.

## Fork

A `fork` is used to run `a set of tasks simultaneously` but **it does not 
wait for the results** from all tasks within `fork`. Instead, it will branch 
the pipeline in the sense that each task within `fork` will have their own 
set of downstream tasks. It

* must **multiply each downstream task** (tasks proceeding the `fork`) as many 
times 
as the `fork` branches (number of tasks within `fork`).
* currently uses `join` operator to create each branch.

*Example*

```javascript
// executes task 1 and task 2 at the same time
// However it does not wait for task1 or task2 to finish before executing 
// task3 after task1 and after task2
const pipeline3 = join(fork(task1, task2), task3)
// executes the pipeline
pipeline3()
```

The above referenced `pipeline3` will run both tasks  (`task1` and `task2`) 
simultaneously and will run downstream task (`task3`) twice, after `task1` and 
after 
`task2`, respectively.

### Conceptual difference between junction and fork

While `junction` is a tree with several branches but that end up in the same 
leaf, `fork` is a tree with several branches that each one end up on its own 
leaf. 

Therefore, `junction` should be used everytime the user wants to 
**wait** for the results from the tasks within `junction` and then perform 
other tasks:

```javascript
const pipeline4 = join(junction(taskA, taskB), taskC)
```

In `pipeline4` example `taskA` and `taskB` will be executed simultaneously but 
then `taskC` will be waiting for the `task.finish` event for both `taskA` and
 `taskB`. This behavior is most useful when `taskC` inputs depend on both 
 tasks (`taskA` and `taskB`) outputs.
 
 On the other hand, `fork` should be used everytime the user **do not want to 
 wait** for all tasks to finish before running downstream tasks:
 
 ```javascript
const pipeline5 = join(fork(task1, task2), task3)
```

In `pipeline5`, if for instance `task1` is finished but `task2` is not, 
`task3` will run in a independent branch after `task1`. Then, when `task2` is
 finished, `taks3` will be run again after `task2`. In this case we will end 
 up with two branches:
 
 ```javascript
taks1 --> task3
task2 --> task3
```

### Concurrency limitations of fork and junction

`junction` and `fork` may start running n number of processes and consume x 
RAM, since they fire several tasks at the same time. Therefore, user must 
have this into account.

Imagine that `fork`(or `junction`) branches the pipeline into 8 independent 
branches:

```javascript
const pipeline6 = join(
  fork(task1, task2, task3, task4, task5, task6, task7, task8), 
  finalTask
)
```

Now imagine that each task uses 1 CPU and you only have 4 CPUs available, 
this will endup consuming more CPUs than we wanted. So, in current API users 
must handle this manually:

```javascript
// first run
const pipeline6 = join(
  fork(task1, task2, task3, task4),
  finalTask
)
// and then
const pipeline6_2 = join(
  fork(task5, task6, task7, task8),
  finalTask
)
```
