# Forkception

Fork is a special case within bionode-watermill because it requires to know 
everything that follows the tasks within the `fork`. In current API there is
 no easy way to get graph shape before actually running the `pipeline.js` 
 scripts. Therefore, `fork` uses `join` to construct independent lineages for
  each task within the `fork`.
  
  Currently `lib/orchestrator/fork.js` just collects the tasks inside `fork` 
  and then everything else is handled by `lib/orchestrator/join.js`, i.e., 
  the branching of the pipeline in several leaves. This branching requires 
  that `fork` *downstream tasks* are multiplied as many times as the number 
  of tasks inside `fork`.
  
  > For instance if fork contains 3 tasks the pipeline will have to branch
  > 3 times from this point on:
  >
  >```javascript
  >const pipeline = join(task0,fork(task1, task2, task3), task4)
  >```
  > This will end up in three different branches (or lineages) with three 
  >leaves:
  > `task0 --> task1 --> task4`, `task0 --> task2 --> task4` and `task0 --> 
  task3 --> task4`.
  >
  > Note: here task4 is a downstream task.
  
  Notice that task4 needed to be triplicated, since it will have different 
  results 
  depending on the results from its parent tasks (`task1`, `task2` or `task3`).
  
  Therefore, if a `fork` is found within a `join` wrapping the all 
  the pipeline a 
  bunch of rules were designed to make fork multiply each *downstream task* 
  as well as be able to properly work all other orchestrators that have their
   own behavior (check [Orchestration](Orchestration.md)).
   
   ## Current limitations
   
   Take the following pipeline as an example:
   
   ```javascript
const pipeline2 = join(
  task0,
  fork(task1, task2), 
  task3, 
  fork(task4, task5), 
  task6)
```
  
  This pipeline represents a further challenge since it requires to 
  duplica the second `fork` itself rather than their contained tasks, as all 
  other rules made for other orchestrators to work inside `fork`. The same is
   true for `junction`, so if in the above example the second `fork` is in fact
    a `junction`, the problem will remain the same.
    
However, there is a workaround this issue that may be used for now:

```javascript
const pipeline2 = join(
  task0,
  fork(
    join(task1, task3, fork(task4, task5)), 
    join(task2, task3, fork(task4, task5))
    ),
  task6
)
```

In the above pice of code I just added `fork(task4, task5)` inside the 
first `fork` in order to be able to get the same structure. To sum up,
the tasks between the two forks (`task3`) as well as the second `fork` (`fork
(task4, task5)`) have to be manually duplicated.
  
  