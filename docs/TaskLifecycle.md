# Task Lifecycle

1. Creating -> new entry in store w/ defaults + input/output/params

2. is resumable -> **on** or **off** -> set status in store -> skip to 7

3. resolve input -> set `task.resolvedInput`
   - from filesystem if first task in pipeline
   - from **collection** otherwise

4. `operation = operationCreator(resolvedInput)` -> { child process, promise, stream } -> set `task.operation`

5. set writable and/or readable of duplexify from operation

6. catch end/finish/close of duplex (end-of-stream)

7. resolve output -> set `task.resolvedOutput`

   - traverse over output, over validators

## 1. Creating

Add a `task` entry to the `tasks` object of the redux store, by applying `Object.assign({}, defaultTask, userTask)` and creating a `uid` based on hashing the input, output, and params. Thus, if an identical `uid` appears, we can escape creating duplicates.

## 2. Check if resumable

Checks `task.resumable`, if **on** skips to step 7, if **off** continue to step 3

## 3. Resolve Input

###  *via c*ollecton

The `collection` is a saved state between tasks. This lets you match a glob pattern to the output of any previous task in the same *task lineage*. Arrays of parameters can introduce new task lineages.

This lets you use files from tasks more than just the previous task.

### Issues to Watch Out For

* file conflicts:
  * unexpected match
  * unexpected multiple matches
    * perhaps enforce only 1 match, unless specified, or throw error otherwise

Should be OK overall if user is considerate, going to write example pipelines to flesh this issues out:

* use same filename (e.g. `reads.sra`) for different species, but because **collection has objects scoped by task lineage parameters**, it only resolves to the correct file
  * should also validate header of resolved file to see if it matches expected specie, etc
* go backwards through state, when conflicts arise (e.g. two tasks produce the same output file name with different contents), throw an error (user needs to write pipeline better), or pick the most recent one?
* iterative processing on same file, improving results each time (the trinity example)

### via filesystem

This is how input is resolved for the first task in a pipeline.

TODO option to fallback to this if resolving from collection fails?

TODO becomes somewhat useless once tasks run in their own folder.

## 4. Create Operation

By this point, input has been resolved to the filesystem of the collection. Then we call the `operationCreatore` with `resolvedProps`, where `resolvedProps` has `input` replaced with `resolvedInput`:

```javascript
operation = operationCreator(resolvedProps)
```

This should create a stream, child process, or promise.

## 5. Set Duplexify Readable and/or Writable from Operation

Wraps the `operation` into a duplexify stream. This is so that streams, child processes, and promises are all handled the same way, and so that you could apply stream operations, like forking (multi-write-stream) to promises.

## 6. End of Stream

Catch the completion of the duplexify stream with end-of-stream, or `.on('close')` if it was a child process. This is so we can start resolving output only once the operation is actually finished.

## 7. Resolve Output

Resolve output glob patterns to filesystem. Run validators over resolved absolute paths. If all validators succeed, can emit a `task.finish` event.
