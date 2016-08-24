# bionode-watermill

[![npm version](https://badge.fury.io/js/bionode-waterwheel.svg)](https://badge.fury.io/js/bionode-waterwheel) [![node](https://img.shields.io/badge/node-v6.x-blue.svg)]() [![Build Status](https://travis-ci.org/bionode/bionode-watermill.svg?branch=master)](https://travis-ci.org/bionode/bionode-waterwheel)  [![codecov.io](https://codecov.io/github/bionode/bionode-watermill/coverage.svg?branch=master)](https://codecov.io/github/bionode/bionode-waterwheel?branch=master)

*Watermill: A Streaming Workflow Engine*

[![NPM](https://nodei.co/npm/bionode-waterwheel.png?downloads=true&stars=true)](https://nodei.co/npm/bionode-waterwheel/)

**NOTE: v0.2 is a prerelease**

See this [draft blog post](https://github.com/thejmazz/jmazz.me/blob/master/content/post/ngs-workflow.md) for an introduction to why this tool was developed, specifically with respect to improving upon ideas introduced by [Snakemake](https://bitbucket.org/snakemake/snakemake/wiki/Home) and [Nextflow](http://www.nextflow.io/) while introducing the notion of *streaming data* via [Stream](https://nodejs.org/api/stream.html).

### Who is this tool for?

Waterwheel is for **programmers** who desire an efficient and easy-to-write
methodology for developing complex and dynamic data pipelines, while handling
parallelization as much as possible. Waterwheel is an npm module, and is
accessible by anyone willing to learn a little JavaScript. This is in contrast
to other tools which develop their own DSL (domain specific language), which is
not useful outside the tool. By leveraging the npm ecosystem and JavaScript on
the client, Waterwheel can be built upon for inclusion on web apis, modern web
applications, as well as native applications through
[Electron](http://electron.atom.io/). Look forward to seeing Galaxy-like
applications backed by a completely configurable Node API.

Waterwheel is for **biologists** who understand it is important to experiment
with sample data, parameter values, and tools. Compared to other workflow
systems, the ease of swapping around parameters and tools is much improved,
allowing you to iteratively compare results and construct more confident
inferences. Consider the ability to construct your own
[Teaser](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-015-0803-1)
for *your data* with a *simple syntax*, and getting utmost performance out of
the box.

### API (v0.4)

See [docs](./docs/Task.md).

### Example Genomic Pipeline

Once you have a Node runtime installed on your system, we can begin writing the
"hello world" pipeline. For this we will illustrate how streams and
parrallezation can be leveraged to improve the performance of our pipeline while
also improving readability, atomicity of tool descriptions, and
self-documentation of the pipeline.

See [variant-calling-simple.js](./examples/vc-simple/variant-calling-simple.js)
