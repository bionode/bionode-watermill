# Beginner Walkthrough

Bionode-watermill can be used to create complex bioinformatics pipelines,
but at its core are simple ideas. The first of these is the idea of a *task*.

## Task

A `task` describes the idea that *some input* will be transformed into *some
output*. First, we will consider input and output to both be a file.
Before jumping into JavaScript, let's see what our first task is!

Our first simple pipeline will convert a file with lowercase characters into one
with upper case characters.  On your UNIX command line (on Windows, see
[try.bionode.io][trybionode] or Docker):

[trybionode]: try.bionode.io

```bash
# setting up "raw data"
echo abcdefghijklmnopqrstuvqxyz > alphabet.lowercase
# double check its there!
cat alphabet.lowercase

# pipeline
cat alphabet.lowercase | tr '[:lower:]' '[:upper:]' > alphabet.uppercase
# see the results!
cat alphabet.uppercase
```

Okay, so how to do the same thing with bionode-watermill? We can create a `task`
for it. A `task` is created by passing `watermill.task` two parameters

- object describing `input`, `output`, `params`, `name` ("props")
- the "operation creator", a function which will receive "props"

For now, let's forget about props and just translate the simple pipeline
above into a `task`:

```javascript
// pipeline.js

// Make sure to "npm install bionode-watermill" first
const { task } = require('bionode-watermill')

const uppercaser = task({
  input: '*.lowercase',
  output: '*.uppercase',
  name: 'Uppercase *.lowercase -> *.uppercase'
}, function(props) {
  const input = props.input
  const output = input.replace(/lowercase$/, 'uppercase')

  return `cat ${input} | tr '[:lower:]' '[:upper:]' > ${output}`
})

uppercaser()
  .then(() => console.log('Task finished!'))
  .catch(err => console.error(err))
```

Then run it:

```bash
node pipeline.js
```

What's going on here? Behind the scenes bionode-watermill will start looking for
files matching the [glob][node-glob] pattern `*.lowercase`. The file you created
earlier will be found, *and its absolute path will be stored*. Then a unique
folder to run the task is created within the `data` folder - this prevents us
from worrying about file overwriting in larger pipelines. Within the folder for
the instance of this task, a symlink to the `*.lowercase` file is created.
Finally, `props`, which holds the absolute path of `*.lowercase` inside
`input` (the same as it was for "props", just now the glob pattern is a valid file
path) is passed to our "operation creator". The operation creator returns a String
(using ES6 template literals). When watermill finds a string returned from a task,
it will run it as a shell child process within that task instance's folder.

The log output will have each line prefixed with a hash: this is the ID of the
instance of that task, and is the name of the folder in `data` which was made
for this task, look inside there and you will find the `*.uppercase` file.

[node-glob]: https://github.com/isaacs/node-glob

There are two ways we can improve the readability of our task:

- assignment destructuring

```javascript
const resolvedProps = { input: '/some/path/to/file.txt' }

// Without assignment destructuring
const input = obj.input
// With assignment destructuring
const { input } = obj
```

- "fat arrow" functions

```javascript
// Without =>
function operationCreator(props) {
    return '...'
}.bind(this) // => will automatically ".bind(this)" (so that this is the same inside and outside the function)
// With =>, can return object directly instead of having a function body
const operationCreator = (props) => '...'
```

With those syntax features, our task looks will look like this:

```javascript
const uppercaser = task({
  input: '*.lowercase',
  output: '*.uppercase',
  name: 'Uppercase *.lowercase -> *.uppercase'
}, ({ input }) =>
    `cat ${input} | tr '[:lower:]' '[:upper:]' > ${input.replace(/lowercase$/, 'uppercase')}`
)
```

Sometimes it is more appropriate to use regular functions, or not to use
assignment destructuring, depending on the task. For example, not using `=>` so
we can declare `output` within a function body. I explain these two language
features because further examples may use them.

How about using Node streams to perform the uppercase transformation? If the
operation creator returns a stream, watermill will understand that. Then we can
replace our child processes with a transform stream using the module
[through2][through2] as a helper to create a transform stream. The task looks like
this:

```javascript
const fs = require('fs')

const { task } = require('bionode-watermill')
const through = require('through2')

const uppercaser = watermill.task({
  input: '*.lowercase',
  output: '*.uppercase',
  name: 'Uppercase *.lowercase -> *.uppercase with Node transform stream'
}, ({ input }) =>
  fs.createReadStream(input)
    .pipe(through(function(chunk, enc, next) {
      // don't use => because we want `this` to come from stream
      this.push(chunk.toString().toUpperCase())
      // tell through we are ready for next chunk, pass no error
      next(null)
    }))
    .pipe(fs.createWriteStream(input.replace(/lowercase$/, 'uppercase')))
)
```

[through2]: https://github.com/rvagg/through2#readme

### Exercise

1. Write a task that downloads SRA for a given species using
   [bionode-ncbi][bionode-ncbi]. *HINT: you can use CLI or Node module;
   watermill needs to be given a flowing stream in the operation creator in
   order to watch it finish*

## Orchestrators

The next core idea of bionode-watermill is that of *orchestrators*. 
Orchestrators are
ways to combine tasks. 

### Join

The simplest of the operators is `join`. `join` takes
any number of tasks as parameters, and will run each of the tasks sequentially.
Where the first task will check the current working folder for files that match
`input`, within a join lineage, downstream tasks will attempt to match their
`input` to an upstream task's `output`.

A simple join pipeline is a download + extract one. For this we can download 
some
reads from the NCBI Sequence Read Archive (SRA) and extract them with
`fastq-dump` from sratools.

```javascript
// Define getSamples as a higher order task so we can create different ones
// for different accessions
const getSamples = (accession) => task({
  params: {
    db: 'sra',
    accession: accession
  },
  output: '**/*.sra',
  dir: process.cwd(), // Set dir to resolve input/output from
  name: `Download SRA ${config.sraAccession}`
}, ({ params }) => ncbi.download(params.db, params.accession).resume() )

const fastqDump = task({
  input: '**/*.sra',
  output: [1, 2].map(n => `*_${n}.fastq.gz`),
  name: 'fastq-dump **/*.sra'
}, ({ input }) => `fastq-dump --split-files --skip-technical --gzip ${input}` )

const pipeline = join(getSamples('2492428'), fastqDump)

pipeline()
```

### Junction

The next operator you might use is `junction`. Junction lets you run a set of
tasks in parallel. If `junction` is used within a `join`, the downstream tasks
will have access to the outputs of each of the parallelized items from `junction`.

We can use `junction` to run our pipeline above over multiple samples and 
make a manifest file with all `.fna.gz` files that were downloaded:

```javascript
// task to create manifest file with all .fna.gz files downloaded
const listFiles = task({
  input: '*.fna.gz',    // this is used to symlink files to this task directory
  output: '*.txt',
}, ({ input}) => `ls *.fna.gz > listFiles.txt`
)

const pipeline = join(
  junction(
    join(getSamples('2492428'), fastqDump),
    join(getSamples('1274026'), fastqDump)
  ),
  listFiles
)
```

> For more details on multiple inputs check this [link](MultipleInput.md).

### Fork

The last operator is `fork`. Fork lets you replace one task in your pipeline
with more than one, and have the downstream tasks duplicated for each task you
pass into fork. This is useful for comparing tools, options within a tool or 
even perform a couple of repetitive tasks for many samples from
 a given pipeline.

For example, we can run use fork to fetch many samples using 
`getSamples` and extract each one `fastqDump` and then uncompressed them all 
using out new `gunzipIt` task:

```javascript
const gunzipIt = task({
  input: '*.fna.gz',
  output: '*.fa',
}, ({ input}) => `gunzip -c ${input} > ${input.split('.')[0]}.fa`
)

const pipeline = join(
  fork(
    join(getSamples('2492428'), fastqDump),
    join(getSamples('1274026'), fastqDump)
  ),
  gunzipIt
)
```

What's important to note is that *we only define one
`gunzipIt` task, yet it is ran on each task in the `fork`*. Also, `gunzipIt` 
will start running once each task within the `fork` gets finished and thus it
 does not have to wait for all tasks within `fork` to finish.

## Pipeline execution

Once you have defined your **tasks** and their **orchestration** make sure 
you call at the end of each `pipeline.js`:

```javascript
pipeline() // makes the assembled pipeline executed
``` 

and then on your terminal:

```bash
node pipeline.js
```

You can also use **google chrome inspector** if you want to debug something:

```bash
node --inspect pipeline.js
```

>**[Only recommended for more advanced users]**
>
>You may even go further by using an environmental variable 
>`REDUX_LOGGER=1`, which basically lets you log more information (*redux* 
>*actions* and *states*) on the workflow within each `bionode-watermill` task:
>
>```bash
>REDUX_LOGGER=1 node --inspect pipeline.js
>```

## Summary

You've seen that bionode-watermill can be used to defined `tasks` that *resolve*
their `input` and `output`. These tasks can be child processes (i.e. a program
on the command line) or just regular JavaScript functions. As well, a brief 
introduction
to the operators `join`, `junction`, and `fork` was made.

These are all the tools necessary to make a bioinformatics pipeline! 

Check out our [tutorial](https://github.com/bionode/bionode-watermill-tutorial)
if you want to get a feel on how it looks like and start assembling your own 
pipelines.

Here are some **example pipelines** 

[Variant calling pipeline](https://github.com/bionode/bionode-watermill/blob/master/examples/variant-calling-filtered/pipeline.js):

```javascript
const pipeline = join(
  junction(
    join(getReference, bwaIndex),
    join(getSamples, fastqDump)
  ),
  trim, mergeTrimEnds,
  decompressReference, // only b/c mpileup did not like fna.gz
  join(
    fork(filterKMC, filterKHMER),
    alignAndSort, samtoolsIndex, mpileupAndCall // 2 instances each of these
  )
)
```

[Mapping with bowtie2 and bwa](https://github.com/bionode/bionode-watermill/tree/master/examples/pipelines/two-mappers)

```javascript
const pipeline = join(
  junction(
    getReference,
    join(getSamples,fastqDump)
  ),
  samtoolsFaidx, gunzipIt, /* since this will be common to both mappers, there
   is no
   need to be executed after fork duplicating the effort. */
  fork(
    join(indexReferenceBwa, bwaMapper),
    join(indexReferenceBowtie2, bowtieMapper)
  ),
  /* here the pipeline will break into two distinct branches, one for each
   mapping approach used. */
  samtoolsView, samtoolsSort, samtoolsIndex, samtoolsDepth
)
```

Now go out there and make something cool!

If you would like to help out, see the
[issues](https://github.com/bionode/bionode-watermill/issues).  There is much to
be done, and anything and everything is appreciated.


[bionode-ncbi]: https://github.com/bionode/bionode-ncbi
