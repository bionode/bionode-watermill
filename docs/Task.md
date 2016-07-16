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

`task` takes two parameters: **props** and **actionCreator**:

```javascript
const myTask = task(props, actionCreator)
```

**props** is an object with the following structure:

```javascript
const props = {
  input: 'foo', // valid input, see below
  output: { file: 'bar.txt' }, // valid output, see below
  name: 'My Task',
  alwaysRun: true // other options, see options
}
```

*input* and *output* **are required**. It is fine to run with no task name, a hashed one will be made for you. However, properly named tasks will help greatly reading pipeline output

**actionCreator** is a function that will be provided with a **resolved props object**. `actionCreator` should return a **stream**, **promise**, **curried callback**,  or **value**. If the action creator does not return a stream, it will be wrapped into a stream internal (e.g. `StreamFromPromise`). An example action creator is

```javascript
function actionCreator(resolvedProps) {
   const fs = require('fs')
   const intoStream = require('into-stream')
   const ws = intoStream(resolvedProps.input).pipe( fs.createWriteStream('bar.txt') )
   return ws
}
```

> *Note*
>
> With assignment destructuring and arrow functions, you can write cleaner action creators:

```javascript
const actionCreator = ({ input }) => intoStream(input).pipe(fs.createWriteStream('bar.txt'))
```

## Input and Output

The `input` and `output` objects can be a **value**, **file**, or **stream**  (or an array, object of these). These are **types**.

### Types

#### value
Pass something in directly, as it is, with **no resolution**. Could be a `boolean`, `string`, `number`, `Array`, `Function`, `Promise`, ...

```javascript
{ value: 'foo' }
```

**TODO**: in v0.3, this will be refactored into just

```javascript
'foo'
```

That is, value is the default type, unless a *file* or *stream* is specified.

#### file

a *glob expression*/*regex* to be **resolved to a full path** before/after task


*as glob*:

```javascript
{ file: '**/*.sra' }
```

*as regex*:

```javascript
{ file: /\.bam$/ }
```



**TODO**: potentially enforce only regex, then can be done with the `file` keyed object and just have

```javascript
/\.bam$/
```

at the expense of a somewhat more verbose syntax.

#### stream

Use this if the task should take a stream as input, output, or both.

*task 1*
```javascript
{ input: { value: '2492428' }, output: { stream: 'stdout' } }
```
*task 2*
```javascript
{ input: { stream: 'stdin' }, output: { file: '*_genomic.fna.gz' } }
```
This will take the `output` stream of the preceding task and pipe it into the `input` of this task.

It is also possible to stream out an existing (or as it is created) file:

```javascript
{ input: { 'stream-file': '*.sam' }, output: { stream: 'stdout' } }
```



**TODO**: potentially assume stream if **input and output are not provided**. Then *through tasks* could be defined like:

```javascript
const filterSpecies = task(
  { name: 'Filter species by taxonomy ', params: { taxonomy: 'plantae' } },
  ({ input, params }) => input.pipe(through(function(chunk, enc, cb) {
    if (chunk.taxonomy === params.taxonomy) {
      this.push(chunk)
    }
    cb()
  }))
)
```

## Resolution

Input resolution is  a **reduction over the task props**. 

### Input Resolution



### Output Resolution

The resolved values can be accessed from `myTask.output()` after the task has emitted a `close` event. (**TODO**: verify this, not only that event).

> **NOTE**
>
> Input and output in these examples are not input and output from task. The input is the original output, and the output is the resolved output.

```javascript
// Atomic item
// Input
{ file: '*.sra' }
// Output
'human.sra'

// Array of items
// Input
[{ file: '*.bam' }, { file: '*.fastq.gz' }]
// Output
['reads.bam', 'human.fastq.gz']

// Can mix files and values
// Input
[{ file: '*.bam' }, 'foo', { file: '*.fastq.gz' }]
// Output
['reads.bam', 'foo', 'human.fastq.gz']

// Object of items, arrays of items,  arrays of objects of items, ...
// PS. playing with idea of file being only a regex in `alignment`
// Input
{
  reads: [1, 2].map( num => ({ file: `*_${num}.fastq.gz`}) ),
  alignment: /\.bam$/
}
// Output
{
  reads: ['human_1.fastq.gz', 'human_2.fastq.gz'],
  alignment: 'reads.bam'
}

```



## Task Lifecycle

1. Check if skippable (resolve output and validate)
2. Resolve input via store and/or fs -> `resolvedProps`
3. `action = actionCreator(resolvedProps)` -> stream, promise, curry cb, value
4. action -> readable and/or writable duplexify (if `output !== stream`, wait until end/finish/close, else, pump stream through)
5. resolve output, validate output, store somewhere

### Check if skippable

- traverse output, run validators (existing, non-null, time, hash, custom, etc)

- if it passes:

  - produced the exact same `task.output()` as if it ran

- else:

  - go to the next step

### Resolve Input

####  *via store*

The store is a saved state between tasks. This lets you use files from tasks more than just the previous task. 

```javascript
// store
['human.sra']
// input
'*.sra'
// resolved
'human.sra'

// Only takes what is needed
// store
['human.sra', 'reference.fastq.gz']
// input
'*.fastq.gz'
// resolved
'reference.fastq.gz'

// Can return an array of matches
// store
['reads_1.fastq.gz', 'reads_2.fastq.gz']
// input
'*.fastq.gz'
// resolved
['reads_1.fastq.gz', 'reads_2.fastq.gz']
// In cases like this, try to be as explicit as possible, using input like
[1, 2].map(num => `*_${num}.fastq.gz`)

// Potentially match keys from store, but that hinders task modularity. Instead, store can be reduced down to its filenames only, and input patterns matched to that
// store
{ ref: 'reference.fastq.gz', reads: 'humans.sra' }
// is effectively: ['reference.fastq.gz', 'humans.sra']
// TODO how did store end up with keys like this?
// input
{ ref: /\.fastq.gz$/, reads: { file: '*.sra' } }
// resolved
{ ref: 'reference.fastq.gz', reads; 'human.sra' }
// TODO potentially match input key `ref` to store key `ref`
// Breaks modularity a bit, but also improves specificity to an extent
```

#### Issues to Watch Out For

* file conflicts:
  * unexpected match
  * unexpected multiple matches
    * perhaps enforce only 1 match, unless specified, or throw error otherwise

Should be OK overall if user is considerate, going to write example pipelines to flesh this issues out:

* use same filename (e.g. `reads.sra`) for different species, but because **store has objects hashed by task lineage parameters**, it only resolves to the correct file
  * should also validate header of resolved file to see if it matches expected specie, etc
* go backwards through state, when conflicts arise (e.g. two tasks produce the same output file name with different contents), throw an error (user needs to write pipeline better), or pick the most recent one?
* iterative processing on same file, improving results each time (the trinity example)

#### via filesystem

TODO, may be removed as an option once each task runs in its own folder. This is same as above, but instead of file paths in the store, uses files from current directory.


### Create Action

TODO. Wraps promise into readable stream that emits one chunk for example.

### Duplexify from Action

TODO. Set readable/writable streams asynchronously, as appropiate.

### Resolve Output

TODO. Similar to examples above for input resolution, but to the filesystem. Can run custom validators (e.g. are `reads_1` and `reads_2` correct siblings?).