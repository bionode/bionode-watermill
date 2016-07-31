# Task

The `task` function is the fundamental unit to build pipelines.

It is provided with the default export of `bionode-waterwheel`:

```javascript
const waterwheel = require('bionode-waterwheel')
const task = waterwheel.task
// Or, with assignment destructuring:
const { task } = waterwheel
```

## API

`task` takes two parameters: **props** and **operationCreator**:

```javascript
const myTask = task(props, operationCreator)
```

**props** is an object with the following structure:

```javascript
const props = {
  input: 'foo.txt', // valid input, see below
  output: 'bar.txt', // valid output, see below
  name: 'My Task',
  alwaysRun: true // force rerun even if output already exists
  // other options, see options
}
```

*input* and *output* are required for tasks that deal with files. If either is not provided, it will be assumed the task is then a *streaming task* - i.e., it is a duplex stream with writable and/or readable portions. Consider: 

```javascript
const throughCapitalize = through(function (chunk, env, next) {
  // through = require('through2') - a helper to create through streams
  // takes chunk, its encoding, and a callback to notify when complete pushing
  // push a chunk to the readable portion of this through stream with
  this.push(chunk.toString().toUpperCase())
  // then call next so that next chunk can be handled
  next()
})

const capitalize = task({
  name: 'Capitalize Through Stream'
}, 
// Here, input is a readable stream that came from the previous task
// Let's return a through stream that takes the input and capitalizes it                 
({ input }) => input.pipe(throughCapitalize) )
```

You could connect `capitalize` to a readable and writable file stream with:

```javascript
const readFile = ({
  input: '*.lowercase',
  name: 'Read from *.lowercase'
}, ({ input }) => {
  const rs = fs.createReadStream(input)
  // Add file information to stream object so we have it later
  rs.inFile = input
})

const writeFile = ({
  output: '*.uppercase',
  name: 'Write to *.uppercase'
}, ({ input }) => fs.createWriteStream(input.inFile.swapExt('uppercase')))

// Can now connect the three:
join(readFile, capitalize, writeFile)
```

Of course, this could be written as one single task. This is somewhat simpler, but the power of splitting up the read, transform, and write portions of a task will become apparent once we can provide multiple sets of parameters to the transform and observe the effect, *without having to manually rewire input and output filenames*. As a single task the above would become:

TODO introduce single task version first

```javascript
const capitalize = task({
  input: '*.lowercase',
  output: '*.uppercase',
  name: 'Capitalize *.lowercase -> *.uppercase'
}, ({ input }) =>
	fs.createReadStream(input)
	.pipe(throughCapitalize)
	.pipe(fs.createWriteStream(input.swapExt('lowercase')))                       
)
```

It is fine to run with no task name, a hashed one will be made for you. However, properly named tasks will help greatly reading pipeline output

**operationCreator** is a function that will be provided with a **resolved props object**. `operationCreator` should return a **stream** or a **promise**. If the operation creator does not return a stream, it will be wrapped into a stream internally (e.g. `StreamFromPromise`). An example operation creator is

```javascript
function operationCreator(resolvedProps) {
   const fs = require('fs')
   const intoStream = require('into-stream')
   const ws = intoStream(resolvedProps.input).pipe( fs.createWriteStream('bar.txt') )
   return ws
}
```

> *Note*
>
> With assignment destructuring and arrow functions, you can write cleaner operation creators:

```javascript
const operationCreator = ({ input }) => intoStream(input).pipe(fs.createWriteStream('bar.txt'))
```

## Input and Output

The `input` and `output` objects can be a **string glob pattern**, or a plain object of them. TODO an array will introduce task forking. The glob will be **resolved to an absolute path** when it is passed to the `operationCreator`.

For example,

```javascript
{ input: '**/*.sra' }
```

will resolve to something like:

```javascript
{ input: '/data/ERR1229296.sra' }
```

And

```javascript
{ 
  input: {
    reference: '*_genomic.fna.gz',
    reads: ['*_1.fastq.gz', '*_2.fastq.gz']
  }
}
```

will resolve to something like:

```javascript
{
  input: {
    reference: '/data/GCA_000988525.2_ASM98852v2_genomic.fna.gz',
    reads: ['/data/ERR1229296_1.fastq.gz', '/data/ERR1229296_2.fastq.gz']
  }
}
```

If `input` is not provided, the `operation` will be a duplexify stream with the writable portion set. Similarly, if `output` is not provided, the `operation` will be a duplexify stream with the readable portion set. If neither is provided, both the readable and writable portions will be set: the task becomes a *through task*.

An example *through task*:

```javascript
const filterSpecies = task({
  name: 'Filter species by taxonomy',
  params: { taxonomy: 'plantae' }
},
({ input, params }) => input.pipe(through.obj(function (chunk, enc, next) {
    if (chunk.taxonomy === params.taxonomy) {
      this.push(chunk)
    }
    next()
  }))
)
```

## Resolution

Input resolution is  a **reduction over the task props**. 

### Input Resolution

Match to filesystem if in first task in pipeline, otherwise glob patterns are matched to the **collection**.

### Output Resolution

The resolved values can be accessed from `myTask.resolvedOutput` after the task has emitted a `task.finish` event.

```javascript
// Atomic item
// Output
'*.sra'
// Resolved Output
'/data/human.sra'

// Array of items
// Output
['*.bam', '*.fastq.gz']
// Resolved Output
['/data/reads.bam', '/data/human.fastq.gz']

// Plain object of items
// Output
{
  reads: [1, 2].map(n => `*_${n}.fastq.gz`),
  alignment: '*.bam'
}
// Resolved Output
{
  reads: ['human_1.fastq.gz', 'human_2.fastq.gz'],
  alignment: 'reads.bam'
}
```



## Task Lifecycle

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

### 1. Creating

Add a `task` entry to the `tasks` object of the redux store, by applying `Object.assign({}, defaultTask, userTask)` and creating a `uid` based on hashing the input, output, and params. Thus, if an identical `uid` appears, we can escape creating duplicates.

### 2. Check if resumable

Checks `task.resumable`, if **on** skips to step 7, if **off** continue to step 3

### 3. Resolve Input

####  *via c*ollecton

The `collection` is a saved state between tasks. This lets you match a glob pattern to the output of any previous task in the same *task lineage*. Arrays of parameters can introduce new task lineages.

This lets you use files from tasks more than just the previous task. 

#### Issues to Watch Out For

* file conflicts:
  * unexpected match
  * unexpected multiple matches
    * perhaps enforce only 1 match, unless specified, or throw error otherwise

Should be OK overall if user is considerate, going to write example pipelines to flesh this issues out:

* use same filename (e.g. `reads.sra`) for different species, but because **collection has objects scoped by task lineage parameters**, it only resolves to the correct file
  * should also validate header of resolved file to see if it matches expected specie, etc
* go backwards through state, when conflicts arise (e.g. two tasks produce the same output file name with different contents), throw an error (user needs to write pipeline better), or pick the most recent one?
* iterative processing on same file, improving results each time (the trinity example)

#### via filesystem

This is how input is resolved for the first task in a pipeline.

TODO option to fallback to this if resolving from collection fails?

TODO becomes somewhat useless once tasks run in their own folder.

### 4. Create Operation

By this point, input has been resolved to the filesystem of the collection. Then we call the `operationCreatore` with `resolvedProps`, where `resolvedProps` has `input` replaced with `resolvedInput`:

```javascript
operation = operationCreator(resolvedProps)
```

This should create a stream, child process, or promise.

### 5. Set Duplexify Readable and/or Writable from Operation

Wraps the `operation` into a duplexify stream. This is so that streams, child processes, and promises are all handled the same way, and so that you could apply stream operations, like forking (multi-write-stream) to promises.

### 6. End of Stream

Catch the completion of the duplexify stream with end-of-stream, or `.on('close')` if it was a child process. This is so we can start resolving output only once the operation is actually finished.

### 7. Resolve Output

Resolve output glob patterns to filesystem. Run validators over resolved absolute paths. If all validators succeed, can emit a `task.finish` event.