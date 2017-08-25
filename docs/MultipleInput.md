# Multiple inputs on tasks

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