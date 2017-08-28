# Internals

This section documents how bionode-watermill functions internally. Read this if you would like to understand the tool better, or would like to contribute *middlewares*.

Watermill and [Gulp](https://github.com/gulpjs/gulp) are alike in some ways:

- both provide an interface to run **tasks** in **series and parallel**
- both maintain a **registry of tasks** - Gulp uses [Undertaker](https://github.com/gulpjs/gulp) and Watermill uses a [Redux](https://github.com/reactjs/redux) store
- both use [async-done](https://github.com/gulpjs/async-done) to catch errors and completion from the task function

Yet also differ in some respects:

- Undertaker uses [bach](https://github.com/gulpjs/bach) for composing async 
functions in serial and parallel, Watermill uses its own serial and parallel 
implementations, which differ in one important way: **objects are passed through between tasks**
- Watermill's `task` is not actually a "task function" (i.e. what is sent into 
async-done), but rather a wrapper around that which has an extensive lifecycle, producing the task function around the middle
- Gulp uses [vinyl-fs](https://github.com/gulpjs/vinyl-fs) to store applied 
transformations to a source stream of files (e.g. `gulp.src('*.lowercase')
.pipe(capitalize()).pipe(gulp.dest('uppercase')) )`) while Watermill does not: 
the data sets are too big to be stored in one Buffer
- In addition to serial and parallel, Watermill provides extra operators like 
fork

Watermill was created to enable the development of data **pipelines with 
modular and swappable components**. As a result, **parameters and options are the first class citizens, while file/folder names are handled automatically**. The streaming nature lets you **mix streams and processes**: e.g. filter search results in a through stream, run a pipeline on each value.








