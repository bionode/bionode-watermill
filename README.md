<p align="center">
  <a href="http://bionode.io">
    <img height="200" width="200" title="bionode" alt="bionode logo" src="https://rawgithub.com/bionode/bionode/master/docs/bionode-logo.min.svg"/>
  </a>
  <br/>
  <a href="http://bionode.io/">bionode.io</a>
</p>

# bionode-watermill

> Bionode-watermill: A (Not Yet Streaming) Workflow Engine

[![npm version](https://badge.fury.io/js/bionode-watermill.svg)](https://badge.fury.io/js/bionode-watermill) 
[![node](https://img.shields.io/badge/node-v6.x-blue.svg)]() 
[![Build Status](https://travis-ci.org/bionode/bionode-watermill.svg?branch=dev)](https://travis-ci.org/bionode/bionode-watermill)  
[![codecov.io](https://codecov.io/github/bionode/bionode-watermill/coverage.svg?branch=master)](https://codecov.io/github/bionode/bionode-watermill?branch=master)
[![Gitter](https://img.shields.io/gitter/room/nwjs/nw.js.svg)](https://gitter.im/bionode/bionode-watermill)

[![NPM](https://nodei.co/npm/bionode-watermill.png?downloads=true&stars=true)](https://nodei.co/npm/bionode-watermill/)

## Table of Contents

* [What is bionode-watermill](#what-is-bionode-watermill)
    * [Main features](#main-features)
    * [Who is this tool for?](#who-is-this-tool-for)
* [Installation](#installation)
* [Documentation](#documentation)
* [Tutorial](#tutorial)
* [Example pipelines](#example-pipelines)
* [Why bionode-watermill?](#why-bionode-watermill)
* [Contributing](#contributing)




## What is bionode-watermill

**Bionode-watermill** is a workflow engine that lets you assemble and run 
bioinformatic pipelines with ease and less overhead. Bionode-watermill 
pipelines are 
essentially node.js scripts in which [tasks](docs/BeginnerWalkthrough.md#task) are the modules that will be 
assembled in the final *pipeline* using [orchestrators](docs/BeginnerWalkthrough.md#orchestrators).

### Main features

* Modularity
* Reusability
* Automated Input/Output handling
* Ability to run programs using Unix shell
* Node.js integration
* [Streamable tasks](docs/Task.md#streamable-tasks-potential) (still not 
implemented - Issue [#79](https://github.com/bionode/bionode-watermill/issues/79))

### Who is this tool for?

Bionode-watermill is for **biologists** who understand it is important to 
experiment with sample data, parameter values, and tools. Compared to other 
workflow systems, the ease of swapping around parameters and tools is much 
improved, allowing you to iteratively compare results and construct more 
confident inferences. Consider the ability to construct your own 
[Teaser](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-015-0803-1) 
for *your data* with a *simple syntax*, and getting utmost performance out of the box.


Bionode-watermill is for **programmers** who desire an efficient and 
easy-to-write methodology for developing complex and dynamic data pipelines, 
while handling parallelization as much as possible. Bionode-watermill is an npm 
module, and is accessible by anyone willing to learn a little JavaScript. This 
is in contrast to other tools which develop their own DSL 
(domain specific language), which is not useful outside the tool. By leveraging 
the npm ecosystem and JavaScript on the client, Bionode-watermill can be built 
upon for inclusion on web apis, modern web applications, as well as native 
applications through [Electron](http://electron.atom.io/). Look forward to 
seeing Galaxy-like applications backed by a completely configurable Node API.


## Installation

Local installation:

```npm install bionode-watermill```

Global installation:

```npm install bionode-watermill -g```

## Documentation

Our documentation is available [here](https://thejmazz.gitbooks.io/bionode-watermill/content/). 
There you may find how to **use** bionode-watermill to construct and **run** 
your 
pipelines. Moreover, you will also find the description of the API to help 
anyone 
willing to **contribute**.


## Tutorial

- [Try bionode-watermill tutorial!](https://github.com/bionode/bionode-watermill-tutorial)

## Example pipelines

- [Toy pipeline with shell/node](https://github.com/bionode/bionode-watermill/blob/master/examples/pipelines/pids/pipeline.js)
- [Simple capitalize task](https://github.com/bionode/bionode-watermill/blob/master/examples/pipelines/capitalize/capitalize.js)
- [Simple SNP calling](https://github.com/bionode/bionode-watermill/blob/master/examples/pipelines/variant-calling-simple/pipeline.js)
- [SNP calling with filtering and fork](https://github.com/bionode/bionode-watermill/blob/master/examples/pipelines/variant-calling-filtered/pipeline.js)
- [Mapping with bowtie2 and bwa (with tutorial)](https://github.com/bionode/bionode-watermill/tree/master/examples/pipelines/two-mappers)

## Why bionode-watermill?

[This blog post](https://jmazz.me/blog/NGS-Workflows)
compares the available tools to deal with NGS workflows, explaining the 
advantages of each one, including **bionode-watermill**.


## Contributing
We welcome all kinds of contributions at all levels of experience, please 
refer to 
the [Issues section](https://github.com/bionode/bionode-watermill/issues). 
Also, you can allways reach us on [gitter](https://gitter.im/bionode/bionode-watermill).

### Feel free to submit your pipeline to us

Just make a PR for us, that adds a pipeline under `./examples/pipelines/`. 
You can check some of the already existing examples [here](examples/pipelines).
