# Task

The `task` function is the fundamental unit to build pipelines.

It is provided with the default export of `bionode-watermill`:

```javascript
const watermill = require('bionode-watermill')
const task = watermill.task
// Or, with assignment destructuring:
const { task } = watermill
```

## API

`task` takes two parameters: **props** and **operationCreator**:

```javascript
const myTask = task(props, operationCreator)
```

### props

**props** is an object with the following structure:

#### Input / output tasks

```javascript
const props = {
  input: '*.txt', // valid input, see below
  output: '*.txt', // valid output, see below
  name: 'My Task',
  params: { output: 'bar.txt' }, //the actual output name that can be passed to 
  // operationCreator or even other params that you may wish to pass to 
  // operationCreator.
  alwaysRun: true // force rerun even if output already exists
  // other options, see options
}
```

*input* and *output* patterns are required for tasks that deal with files. 
*params.output* allows to name the output files, while *output* is used to 
check if output was properly resolved.

```javascript
// example task with input/output files
const task = ({
  input: '*.txt',
  output: '*.txt',
  params: {
    output: 'bar.txt'
  }
}, ({ input, params }) => `cp ${input} ${params.output}`
)
```

## Streamable tasks

If either (input or output) is 
not provided, it will be assumed the task is then a *streaming task* - i.e., it 
is a duplex stream with writable and/or readable portions. Consider this 
javascript function:

```javascript
const throughCapitalize = through(function (chunk, env, next) {
  // through = require('through2') - a helper to create through streams
  // takes chunk, its encoding, and a callback to notify when complete pushing
  // push a chunk to the readable portion of this through stream with
  this.push(chunk.toString().toUpperCase())
  // then call next so that next chunk can be handled
  next()
})
```

And the following tasks:

```javascript
const capitalize = task({
  name: 'Capitalize Through Stream'
},
// Here, input is a readable stream that came from the previous task
// Let's return a through stream that takes the input and capitalizes it
({ input }) => input.pipe(throughCapitalize) )

const readFile = task({
  input: '*.lowercase',
  name: 'Read from *.lowercase'
}, ({ input }) => {
  const rs = fs.createReadStream(input)
  // Add file information to stream object so we have it later
  rs.inFile = input
})

const writeFile = task({
  output: '*.uppercase',
  name: 'Write to *.uppercase'
}, ({ input }) => fs.createWriteStream(input.inFile.swapExt('uppercase')))
```

You could connect `capitalize` to a readable (`readFile`) and writable 
(`writeFile`) file 
stream with:

```javascript
// Can now connect the three:
join(readFile, capitalize, writeFile)
```

Of course, this could be written as one single task. This is somewhat simpler, 
but the power of splitting up the read, transform, and write portions of a task 
will become apparent once we can provide multiple sets of parameters to the 
transform and observe the effect, *without having to manually rewire input and 
output filenames*. As a single task the above would become:

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

### operationCreator

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
>
>```javascript
>const operationCreator = ({ input }) => intoStream(input).pipe(fs.createWriteStream('bar.txt'))
>```

### Input and Output

The `input` and `output` objects can be a **string glob pattern**, or a plain 
object of them. The glob will be 
**resolved to an absolute path** when it is passed to the `operationCreator`.

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
    reference: '/data/<uid>/GCA_000988525.2_ASM98852v2_genomic.fna.gz',
    reads: ['/data/<uid>/ERR1229296_1.fastq.gz', '/data/ERR1229296_2.fastq.gz']
  }
}
```

If `input` is not provided, the `operation` will be a duplexify stream with the 
writable portion set. Similarly, if `output` is not provided, the `operation` 
will be a duplexify stream with the readable portion set. If neither is 
provided, both the readable and writable portions will be set: the task becomes 
a *through task*.

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

### Resolution

Input resolution is  a **reduction over the task props**.

#### Input Resolution

Match to filesystem if in first task in pipeline, otherwise glob patterns are matched to the **collection**.

#### Output Resolution

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
