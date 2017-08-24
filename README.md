# bionode-watermill

[![npm version](https://badge.fury.io/js/bionode-watermill.svg)](https://badge.fury.io/js/bionode-watermill) [![node](https://img.shields.io/badge/node-v6.x-blue.svg)]() [![Build Status](https://travis-ci.org/bionode/bionode-watermill.svg?branch=dev)](https://travis-ci.org/bionode/bionode-watermill)  [![codecov.io](https://codecov.io/github/bionode/bionode-watermill/coverage.svg?branch=master)](https://codecov.io/github/bionode/bionode-watermill?branch=master)

*Watermill: A Streaming Workflow Engine*

[![NPM](https://nodei.co/npm/bionode-watermill.png?downloads=true&stars=true)](https://nodei.co/npm/bionode-watermill/)


- [CWL?](#cwl)
- [What is a task?](#what-is-a-task)
- [What are orchestrators?](#what-are-orchestrators)
- [Check out bionode-watermill tutorial!](#check-out-bionode-watermill-tutorial)
- [Example pipelines](#example-pipelines)
- [Why bionode-watermill?](#why-bionode-watermill)
- [Who is this tool for?](#who-is-this-tool-for)

Watermill lets you *orchestrate* **tasks** using operators like **join**, **junction**, and **fork**. Each task has a [lifecycle](https://thejmazz.gitbooks.io/bionode-watermill/content/TaskLifecycle.html) where

1. Input [glob patterns](https://github.com/isaacs/node-glob) are resolved to absolute file paths (e.g. `*.bam` to `reads.bam`)
2. The **operation** is ran, passed resolved input, params, and other props
3. The operation completes.
4. Output glob patterns are resolved to absolute file paths.
5. Validators are ran over the output. Check for non-null files, can pass in custom validators.
6. Post-validations are ran. Add task and output to DAG.

## CWL?

Coming soon.

## What is a task?

A `task` is the fundamental unit pipelines are built with. For more details, see [Task](https://thejmazz.gitbooks.io/bionode-watermill/content/Task.html). At a glance, a task is created by passing in **props** and an **operationCreator**, which will later be called with the resolved input. Consider this task which takes a "lowercase" file and creates an "uppercase" one:

```javascript
const uppercase = task({
  input: '*.lowercase',
  output: '*.uppercase'
}, function(resolvedProps) {
  const input = resolvedProps.input

  return fs.createReadStream(input)
  	.pipe(through(function(chunk, enc, next) {
      next(null, chunk.toString().toUpperCase())
  	})
    .pipe(fs.createWriteStream(input.replace(/lowercase$/, 'uppercase')))
})
```

A "task declaration" like above will not immediately run the task. Instead, the task declaration returns an "invocable task" that can either be called directly or used with an orchestration operator. Tasks can also be created to **run shell programs**:

```javascript
const fastqDump = task({
  input: '**/*.sra',
  output: [1, 2].map(n => `*_${n}.fastq.gz`),
  name: 'fastq-dump **/*.sra'
}, ({ input }) => `fastq-dump --split-files --skip-technical --gzip ${input}` )
```

## What are orchestrators?

Orchestrators are functions which can take tasks as params in order to let you compose your pipeline from a high level view. This **separates task order from task declaration**. For more details, see [Orchestration](https://thejmazz.gitbooks.io/bionode-watermill/content/Orchestration.html). At a glance, here is a complex usage of `join`, `junction`, and `fork`:

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

## Check out bionode-watermill tutorial!

- [Try out bionode-watermill tutorial](https://github.com/bionode/bionode-watermill-tutorial)

## Example pipelines

- [Toy pipeline with shell/node](https://github.com/bionode/bionode-watermill/blob/master/examples/pipelines/pids/pipeline.js)
- [Simple capitalize task](https://github.com/bionode/bionode-watermill/blob/master/examples/pipelines/capitalize/capitalize.js)
- [Simple SNP calling](https://github.com/bionode/bionode-watermill/blob/master/examples/pipelines/variant-calling-simple/pipeline.js)
- [SNP calling with filtering and fork](https://github.com/bionode/bionode-watermill/blob/master/examples/pipelines/variant-calling-filtered/pipeline.js)
- [Mapping with bowtie2 and bwa](https://github.com/bionode/bionode-watermill/tree/master/examples/pipelines/two-mappers)

## Why bionode-watermill?

[This blog post](https://jmazz.me/blog/NGS-Workflows)
compares the available tools to deal with NGS workflows, explaining the 
advantages of each one, including **bionode-watermill**.

## Who is this tool for?

Bionode-watermill is for **programmers** who desire an efficient and easy-to-write methodology for developing complex and dynamic data pipelines, while handling parallelization as much as possible. Bionode-watermill is an npm module, and is accessible by anyone willing to learn a little JavaScript. This is in contrast to other tools which develop their own DSL (domain specific language), which is not useful outside the tool. By leveraging the npm ecosystem and JavaScript on the client, Bionode-watermill can be built upon for inclusion on web apis, modern web applications, as well as native applications through [Electron](http://electron.atom.io/). Look forward to seeing Galaxy-like applications backed by a completely configurable Node API.

Bionode-watermill is for **biologists** who understand it is important to experiment with sample data, parameter values, and tools. Compared to other workflow systems, the ease of swapping around parameters and tools is much improved, allowing you to iteratively compare results and construct more confident inferences. Consider the ability to construct your own [Teaser](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-015-0803-1) for *your data* with a *simple syntax*, and getting utmost performance out of the box.
