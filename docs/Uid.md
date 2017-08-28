# Uid

In bionode-watermill we use **uids**, that are basically hashes. These `uids`
 are used to control if a task is already run or not. For instance, imagine 
 you have ran part of your pipeline before and it failed somewhere in the 
 middle. When you attempt to run your `pipeline.js` script again, 
 bionode-watermill will know that some of your tasks have already run and 
 will not run them again. Instead it will be able to crawl into the `./data` 
 folder too look for the inputs required for the tasks that still need to be 
 run.
 
 ## uid at the task level
 
 Uid at the task level is generate in two different ways:
* First, it generated given the `props` (`input`, `output` and `params`) 
passed by 
the user to the task 
definition (this is generated under `lib/reducers/task.js`). Let's call this 
`uid` the *task default* `uid` (before running).

* Second, if we want a task to be ran twice in the same pipeline it cannot 
have the same `uid` otherwise bionode-watermill will not allow the second 
execution of the pipeline. However, it can properly solve this because there 
is a second level of `uid` generation while the task is running. Therefore, 
a task `uid` is modified on running to get its *final or run* `uid`. This new
 `uid` is generated taking into account *task default* `uid` and its parent 
 tasks `uids` and making a unique hash of all these `uids`. This renders that
  the same task can be ran twice in the pipeline **if** their `trajectory` is
   different, i.e., **if** they have different parent tasks.
   
## uid at the orchestrator level

Orchestrators also have uids that are basically an hash of the `uids` of the 
tasks contained in them. These `uids`are used in `lib/orchestrators/join.js` to 
control the `uid` generation for `fork` because each *downstream task* after a 
`fork` must be multiplied as many times as the tasks inside fork (for further 
details on this see [Forkception](Forkception.md)).
