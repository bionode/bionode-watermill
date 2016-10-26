# Beginner Walkthrough

bionode-watermill can be used to create complicated bioinformatics pipelines,
but at its core are simple ideas. The first of these is the idea of a *task*.

A `task` describes the idea that *some input* will be transformed into *some
output*. For our cases, we will consider input and output to both be a file.
Before jumping into the JavaScript, let's see what our first task is!

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

Okay, so how to the same thing with bionode-watermill? We can create a `task`
for it. A `task` is created by passing `watermill.task` two parameters

- object describing `input`, `output`, `name`
- the "operation creator", a function which will receive "resolved props"

For now, let's forget about resolved props and just translate the simple pipeline
above into a `task`:

```javascript
// pipeline.js

// Make sure to "npm install bionode-watermill" first
const watermill = require('bionode-watermill')

const uppercaser = watermill.task({
  input: '*.lowercase',
  output: '*.uppercase',
  name: 'Uppercase *.lowercase -> *.uppercase'
}, function(resolvedProps) {
  const input = resolvedProps.input
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
files matching the [glob][node-glob] pattern `*.alphabet`. The file you created
earlier will be found, *and its absolute path will be stored*. Then a unique
folder to run the task is created within the `data` folder - this prevents us
from worrying about file overwriting in larger pipelines. Within the folder for
the instance of this task, a symlink to the `*.alphabet` file is created.
Finally, `resolvedProps`, which holds the absolute path of `*.alphabet` inside
`input` (the same as it was for "props", just now the glob pattern is a valid file
path) is passed to our "operation creator". The operation creator returns a String
(using ES6 template literals). When watermill finds a string returned from a task,
it will run it as a shell child process within that task instance's folder.

The log output will have each line prefixed with a hash: this is the ID of the
instance of that task, and is the name of the folder in `data` which was made
for this task, look inside there and you will find the `alphabet.uppercase` file.

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
function actionCreator(resolvedProps) {
    return '...'
}.bind(this) // => will automatically ".bind(this)" (so that this is the same inside and outside the function)
// With =>, can return object directly instead of having a function body
const actionCreator = (resolvedProps) => '...'
```

With those syntax features, our task looks like:

```javascript
const uppercaser = watermill.task({
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
    })
    .pipe(fs.createWriteStream(input.replace(/lowercase$/, 'uppercase')))
)
```

[through2]: https://github.com/rvagg/through2#readme

### Exercise

1. Write a task that downloads SRA for a given species using
   [bionode-ncbi][bionode-ncbi]. *HINT: you can use CLI or Node module; watermill
   needs to be given a flowing stream in the operation creator in order to watch it finish*


[bionode-ncbi]: https://github.com/bionode/bionode-ncbi
