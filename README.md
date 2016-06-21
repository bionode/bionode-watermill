# bionode-waterwheel
Waterwheel: A Streaming Workflow Engine

See this [draft blog post](https://github.com/thejmazz/jmazz.me/blob/master/content/post/ngs-workflow.md) for an introduction to why this tool was developed, specifically with respect to improving upon ideas introduced by [Snakemake](https://bitbucket.org/snakemake/snakemake/wiki/Home) and [Nextflow](http://www.nextflow.io/) while introducing the notion of *streaming data* via [Stream](https://nodejs.org/api/stream.html).

### Who is this tool for?

Waterwheel is for **programmers** who desire an efficient and easy-to-write methodology for developing complex and dynamic data pipelines, while handling parrelization as much as possible. Waterwheel is an npm module, and is accessible by anyone willing to learn a little JavaScript. This is in contrast to other tools which develop their own DSL (domain specific language), which is not useful outside the tool. By leveraging the npm ecosystem and JavaScript on the client, Waterwheel can be built upon for inclusion on web apis, and modern web applications, as well as native applications through [Electron](http://electron.atom.io/).

Waterwheel is for **biologists** who understand it is important to experiment with sample data, parameter values, and tools. Compared to other tools, the ease of swapping around paramaters, and tools is much improved, allowing you to iteratively compare results and construct more confident inferences. Consider the ability to construct your own [Teaser](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-015-0803-1) for *your data* with a *simple syntax*, and getting utmost performance out of the box.

### How does it work?

Waterwheel is a pipeline tool following the dataflow paradigm. Data is transferred between tasks using Node streams. This lets us describe UNIX pipes with multiple *tasks* that can then be piped into each other.

### Getting Started

Once you have a Node runtime installed on your system, we can begin writing the "hello world" pipeline. For this we will illustrate how streams and parrallezation can be leveraged to improve the performance of our pipeline while also improving readability, atomicity of tool descriptions, and self-documentation of the pipeline.

```javascript
const waterwheel = require('bionode-waterwheel')
const { shell } = require('bionode-wrapper')
const ncbi = require('bionode-ncbi')

const { task, run } = waterwheel
const { stdout, stdin } = waterwheel.types

// As an sra file is written, it will be piped as stdout to the next task
const samples = task({
  input: {
    accession: 'SRP061902'
  },
  output: stdout('**/*.sra')
}, (props) => ncbi.download.sra(props.input.accession))

// Takes a streaming sra and outputs fastq
const dump = task({
  input: stdin(),
  output: stdout([1, 2].map(num => `*_${num}.fastq.gz`))
}, (props) => shell('fastq-dump --split-files -'))

run(samples, dump)
```
