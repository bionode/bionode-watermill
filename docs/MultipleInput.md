# Multiple inputs

There are two ways multiple inputs can be given to tasks - within the **task 
scope** or within the **pipeline scope**.

## Within task scope

Within task scope multiple input handling can take into account 
bionde-watermill capacities to fetch files using glob patterns (e.g. `*.fas`).
For instance, consider the following task:

```javascript
const concatTask = task({
  input: '*.fas',
  output: '*.fas',
  params: {output: 'concat.fas'}
}, ( object ) => {
  console.log('input_directory: ', object.dir)
  console.log('input_variable: ', object, object.params)
  return `cat ${object.input.join(' ')} > ${object.params.output}`
  }
)
```

This `concatTask` has the ability to fetch several `*.fas` that are present 
in the current working directory and pass it to the `operationCreator` 
function in task (for further reference refer to [Task](Task.md)). However, 
users need to be careful because these glob patterns are not the same as 
shell wildcards. So they must be handled by javascript before passing them to
 the shell:
  
  `object.input.join(' ')` - will transform the array of inputs into a string
   where each `.fas` file will be separated by a space.
   
## Within the pipeline scope

Users may also have multiple samples that want to run through all the 
pipeline. Imagine you have `n` fasta files and you want all of them to go 
through the pipeline. Currently this will require some javascript around 
tasks and pipelines.
See this example:

First lets define a function that executes the desired pipeline for each 
input file or sample:

```javascript
// wrap every task and pipeline of a file in a single function (in this case
// 'pipeline' function
const fileFlow = (infile) => {
  const concatTask = task({
    // notice that the input file (infile) is here passed as a param of the
    // pipeline function
      input: infile,
      output: '*.fas',
      params: {output: 'concat.fas'}
    }, ( object ) => {
      console.log('input_directory: ', object.dir)
      console.log('input_variable: ', object, object.params)
      return `cat ${object.input} > ${object.params.output}`
    }
  )

  const task0 = task({name: 'coco'}, () => `echo "something0"`)

  const task1 = task({name: 'nut'}, () => `echo "something1"`)

  const task2 = task({name: 'foo'}, () => `echo "something2"`)

  // this returns the pipeline for a file
  const pipelineAssemblage = join(concatTask, task0, fork(task1, task2))
  return pipelineAssemblage
}
```

Notice that it basically wrap task definition and pipeline call 
(`pipelineAssemblage`) and returns it.

Then is just a matter of cycling through each file in current working 
directory (for example):

```javascript
// checks for files in current working directory (where the script is executed)
fs.readdir(process.cwd(), (err, files) => {
  files.forEach(infile => {
    // checks file extension, in this case '.fas'
    if (path.extname(infile) === '.fas') {
      // executes the pipeline function for each infile
      const pipelineMaster = fileFlow(infile)
      pipelineMaster()
    }
  })
})
```

This will execute the function that wraps the tasks and pipeline 
(`pipelineAssemblage`) and executes it each time a input file is given to 
this function. This will render one graph for each input file provided for 
the pipeline function (`fileFlow`).

### Concurrency 

As you may imagine such behavior without any control may end up using more 
resources (either CPU or memory) than we have. So one quick way to solve 
this, would be to use `bluebird`:

```javascript
const Promise = require('bluebird')
``` 

You can even pass a argument to be read from the shell:
```javascript
// makes concurrency the third argument of code in shell
const concurrency = parseFloat(process.argv[2] || "Infinity")
```

And call it like this:

```shell
node pipeline.js 2
```

This will tell that you want to have 2 files running at the same time.
Then is just a matter of returning a `bluebird` `Promise.map` with 
`concurrency`:

```javascript
// checks for files in current working directory (where the script is executed)
fs.readdir(process.cwd(), (err, files) => {
  const desiredFiles = []
  files.forEach(infile => {
    // checks file extension, in this case '.fas' and creates an array with
    // just these files
    if (path.extname(infile) === '.fas') {
      desiredFiles.push(infile)
    }
  })
  // bluebird Promise.map is used here for concurrency where we can
  // limit the flow of files into the pipeline. So, if concurrency is set
  // to 1 it will first execute the first file, then the second and so on.
  // But if 2 it will execute the first 2 files and then the other 2.
  return Promise.map(desiredFiles, (infile) => {
    const pipelineMaster = fileFlow(infile)
    return pipelineMaster()
  }, {concurrency})
})
```

Note that it is quite similar to the example above but returns a `Promise
.map`, which also return the pipeline executions, instead of just executing the 
pipeline for each input file. This is one simple way to empower the user to 
control the resources used by their scripts.

### Pratical example

Imagine we have a couple of samples we want to get from NCBI SRA and map them
 against a given reference using both `bowtie2` and `bwa`.
 
 So first we have to define a task that needs to fetch the reference data and
  this task must be run just once:
  
  ```javascript
  // some config variables that will be used in this pipeline
const config = {
  name: 'Streptococcus pneumoniae',
  sraAccession: ['ERR045788', 'ERR016633'],
  referenceURL: 'http://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/007/045/GCF_000007045.1_ASM704v1/GCF_000007045.1_ASM704v1_genomic.fna.gz'
}

// === TASKS ===

// first lets get the reference genome for our mapping
const getReference = task({
  params: { url: config.referenceURL },
  output: '*_genomic.fna.gz',
  dir: process.cwd(),   // this sets directory for outputs!
  name: `Download reference genome for ${config.name}`
}, ({ params, dir }) => {
  const { url } = params
  const outfile = url.split('/').pop()

  // essentially curl -O
  return request(url).pipe(fs.createWriteStream(outfile))
})
```

Then what we need is a way for all other tasks to occur after `getReference` 
task:

```javascript
// task that defines how we get samples
//then get samples to work with
const getSamples = (sraAccession) => task({
  params: {
    db: 'sra',
    accession: sraAccession
  },
  output: '**/*.sra',
  dir: process.cwd(), // Set dir to resolve input/output from
  name: `Download SRA ${sraAccession}`
}, ({ params }) => `bionode-ncbi download ${params.db} ${params.accession}`
)

// ...other tasks...
// the pipeline can be then called like this
getReference().then(results => {
  const pipeline = (sraAccession) => join(
    getSamples(sraAccession),
    fastqDump,
    gunzipIt,
    fork(
      join(indexReferenceBwa, bwaMapper),
      join(indexReferenceBowtie2, bowtieMapper)
    )
  )
// then fetches the samples and executes the remaining pipeline
  for (const sra of config.sraAccession) {
    const pipelineMaster = pipeline(sra)
    pipelineMaster().then(results => console.log("Results: ", results))
  }
})
```

Notice how the pipeline is ran twice (given that we have two inputs in an 
array (`config.sraAccession`)). And notice that we have passed an argument to
`getSamples` task which is each `config.sraAccession`.

This behavior will result in one vertex with the `getReference` task, which 
outputs become available to  each `pipeline` that each sra sample triggers. 
So, we end up with a `pipeline` for each sample as well (two in this case).

> This is possible because we rely on the `getReference` task to save the 
outputs on current woring directory. Then the api searches first for outputs 
that are generated within each bionode-watermill pipeline `currentCollection`
 and then in the current working directory for a matching file (with the 
 expected glob pattern).
>
> This matching is done by the `matchee` function in 
`lib/lifecycle/resolve-input.js`.
>
> The fact that bionode-watermill also searches for inputs under current 
working directory is what allows to run a pipeline after another and still 
get reference to the previous pipeline (in the example given above). 
Therefore, use it as you please but take into account that **bionode-watermill 
cannot map every folder within your file system** and **`currentCollection` 
just 
saves reference to the files that were run within a pipeline**.
