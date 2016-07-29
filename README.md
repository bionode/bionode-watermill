# bionode-waterwheel

[![npm version](https://badge.fury.io/js/bionode-waterwheel.svg)](https://badge.fury.io/js/bionode-waterwheel) [![Build Status](https://travis-ci.org/bionode/bionode-waterwheel.svg?branch=master)](https://travis-ci.org/bionode/bionode-waterwheel) [![codecov.io](https://codecov.io/github/bionode/bionode-waterwheel/coverage.svg?branch=master)](https://codecov.io/github/bionode/bionode-waterwheel?branch=master) 

*Waterwheel: A Streaming Workflow Engine*

**NOTE: v0.2 is a prerelease**

See this [draft blog post](https://github.com/thejmazz/jmazz.me/blob/master/content/post/ngs-workflow.md) for an introduction to why this tool was developed, specifically with respect to improving upon ideas introduced by [Snakemake](https://bitbucket.org/snakemake/snakemake/wiki/Home) and [Nextflow](http://www.nextflow.io/) while introducing the notion of *streaming data* via [Stream](https://nodejs.org/api/stream.html).

### Who is this tool for?

Waterwheel is for **programmers** who desire an efficient and easy-to-write methodology for developing complex and dynamic data pipelines, while handling parallelization as much as possible. Waterwheel is an npm module, and is accessible by anyone willing to learn a little JavaScript. This is in contrast to other tools which develop their own DSL (domain specific language), which is not useful outside the tool. By leveraging the npm ecosystem and JavaScript on the client, Waterwheel can be built upon for inclusion on web apis, modern web applications, as well as native applications through [Electron](http://electron.atom.io/). Look forward to seeing Galaxy-like applications backed by a completely configurable Node API.

Waterwheel is for **biologists** who understand it is important to experiment with sample data, parameter values, and tools. Compared to other workflow systems, the ease of swapping around parameters and tools is much improved, allowing you to iteratively compare results and construct more confident inferences. Consider the ability to construct your own [Teaser](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-015-0803-1) for *your data* with a *simple syntax*, and getting utmost performance out of the box.

### API (v0.2)

See [docs](./docs/Task.md).

### How does it work?

Waterwheel is a pipeline tool following the dataflow paradigm. Data is transferred between tasks using Node streams. This lets us describe UNIX pipes with multiple *tasks* that can then be piped into each other. Why is this is useful? Consider the following
Nextflow process:

```groovy
process alignAndSort {
  container: 'needsManyTools'

  input:
    file reference from referenceGenome
    file referenceIndex from referenceIndex
    file sample from reads
  output: file 'reads.bam' into readsBAM

  """
  bwa mem -t $THREADS $reference $sample | samtools view -Sbh - | samtools sort $sam -o reads.bam
  """
}
```

While this is performant and neat, and in practice most workflows will execute these three commands in this order, there a few issues:
- will need to make a custom Dockerfile, whereas Biodocker has existing images for the individual tools used
- paramaters for multiple tools are bundled into one task. this decreases modularity and experimentation.

The Snakemake counterpart would have the same issue.

Wouldn't it be neat to declare the parameters for `bwa mem` elsewhere than you do for `samtools sort`? With Node we can use the Streams API to make this possible:

```js
const bwaMem = task({
  input: {
    reference: { file: '*_genomic.fna.gz' },
    reads: { file: '*.fastq.gz' }
  },
  output: 'stdout'
}, ({ input }) => shell(`bwa mem ${input.reference} ${input.reads}`))

const samtoolsView = task({
  input: 'stdin',
  output: 'stdout'
}, () => shell(`samtools view -Sbh -`))

const samtoolsSort = task({
  input: 'stdin'
  output: { file: 'reads.bam' }
}, () => shell('samtools sort - -o reads.bam'))

const alignAndSort = join(bwaMem, samtoolsView, samtoolsSort)
```

This way, we still maintain the "complete task", yet each component can be
configured  on its own. With a public repository of tools, this would enable
increased modularity and the ability to plug and play together pipelines.
Moreover, each task can be  easily Dockerized without the use of a homemade, one
time use, polyglot Dockerfile.

`bwa mem` has many options. None are actually being used above. But consider
something like this:

```js
const bwaMems = ['-option1', '-option2'].map(opts => makeBwaMemTask(opts))

const alignAndSortAndCompareOptions = join(bwaMems, samtoolsView, samtoolsSort)
```

In this case, `bwaMems` will return multiple tasks. Join will automatically
handle that and fork the stream upstream to account for this, and create two
(samtoolsView, samtoolsSort) throughways. This will let researchers easily swap
in and out new parameters and tools so they can quickly compare results: tools
act differently on different datasets!

*While `bwa mem | samtools view | samtools sort` may not be the best example for
this, the idea is that one can modularize a whole pipeline down to individual
processes with individual parameters. Then let Waterwheel compose it all together.*

### Example Genomic Pipeline

Once you have a Node runtime installed on your system, we can begin writing the "hello world" pipeline. For this we will illustrate how streams and parrallezation can be leveraged to improve the performance of our pipeline while also improving readability, atomicity of tool descriptions, and self-documentation of the pipeline.

See [variant-calling-simple.js](./examples/vc-simple/variant-calling-simple.js)
