# bionode-waterwheel

[![Build Status](https://travis-ci.org/bionode/bionode-waterwheel.svg?branch=master)](https://travis-ci.org/bionode/bionode-waterwheel) [![codecov.io](https://codecov.io/github/bionode/bionode-waterwheel/coverage.svg?branch=master)](https://codecov.io/github/bionode/bionode-waterwheel?branch=master)

*Waterwheel: A Streaming Workflow Engine*

See this [draft blog post](https://github.com/thejmazz/jmazz.me/blob/master/content/post/ngs-workflow.md) for an introduction to why this tool was developed, specifically with respect to improving upon ideas introduced by [Snakemake](https://bitbucket.org/snakemake/snakemake/wiki/Home) and [Nextflow](http://www.nextflow.io/) while introducing the notion of *streaming data* via [Stream](https://nodejs.org/api/stream.html).

### Who is this tool for?

Waterwheel is for **programmers** who desire an efficient and easy-to-write methodology for developing complex and dynamic data pipelines, while handling parallelization as much as possible. Waterwheel is an npm module, and is accessible by anyone willing to learn a little JavaScript. This is in contrast to other tools which develop their own DSL (domain specific language), which is not useful outside the tool. By leveraging the npm ecosystem and JavaScript on the client, Waterwheel can be built upon for inclusion on web apis, modern web applications, as well as native applications through [Electron](http://electron.atom.io/). Look forward to seeing Galaxy-like applications backed by a completely configurable Node API.

Waterwheel is for **biologists** who understand it is important to experiment with sample data, parameter values, and tools. Compared to other workflow systems, the ease of swapping around parameters and tools is much improved, allowing you to iteratively compare results and construct more confident inferences. Consider the ability to construct your own [Teaser](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-015-0803-1) for *your data* with a *simple syntax*, and getting utmost performance out of the box.

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
    reference: new File('*_genomic.fna.gz'),
    reads: new File('*.fastq.gz')
  },
  output: stdout()
}, ({ input }) => shell(`bwa mem ${input.reference} ${input.reads}`))

const samtoolsView = task({
  input: stdin(),
  output: stdout()
}, () => shell(`samtools view -Sbh -`))

const samtoolsSort = task({
  input: stdin()
  output: new File('reads.bam')
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

### Getting Started

Once you have a Node runtime installed on your system, we can begin writing the "hello world" pipeline. For this we will illustrate how streams and parrallezation can be leveraged to improve the performance of our pipeline while also improving readability, atomicity of tool descriptions, and self-documentation of the pipeline.

```javascript
const waterwheel = require('bionode-waterwheel')
const { task, taskYAML, join, parallel } = waterwheel
const { File } = waterwheel.types
const { shell, shellPipe } = waterwheel.wrappers

const fs = require('fs')
const ncbi = require('bionode-ncbi')
const request = require('request')

const THREADS = 4
const config = {
  sraAccession: '2492428',
  referenceURL: 'http://ftp.ncbi.nlm.nih.gov/genomes/all/GCA_000988525.2_ASM98852v2/GCA_000988525.2_ASM98852v2_genomic.fna.gz'
}

// task(props, cb)
// cb(props) will be called after Files in props are resolved with globby
// cb is the "task action creator" and should return a Stream, Promise, or callback

const samples = task({
  input: {
    db: 'sra',
    accession: config.sraAccession
  },
  output: '**/*.sra'
}, ({ input }) => ncbi.download(input.db, input.accession) )

const fastqDump = task({
  input: new File('**/*.sra'),
  output: [1, 2].map(n => new File(`*_${n}.fastq.gz`))
}, ({ input }) => shell(`fastq-dump --split-files --skip-technical --gzip ${input}`) )

const downloadReference = task({
  input: config.referenceURL,
  output: new File(config.referenceURL.split('/').pop())
}, ({ input, output }) => request(input).pipe(fs.createWriteStream(output.value)) )

const bwaIndex = task({
  input: new File('*_genomic.fna.gz'),
  // Self-document by clearly describing expected output files
  output: ['amb', 'ann', 'bwt', 'pac', 'sa'].map(suffix => new File(`*_genomic.fna.gz.${suffix}`))
}, ({ input }) => shell(`bwa index ${input}`) )

const alignAndSort = task({
  input: [new File('*_genomic.fna.gz'), new File('*.fastq.gz')],
  output: new File('reads.bam')
}, ({ input }) => shellPipe([
  `bwa mem -t ${THREADS} ${input[0]} ${input[1]} ${input[2]}`,
  'samtools view -Sbh -',
  'samtools sort - -o reads.bam'
]) )

const samtoolsIndex = task({
  input: new File('*.bam'),
  output: new File('*.bam.bai')
}, ({ input }) => shell(`samtools index ${input}`) )

const decompressReference = task({
  input: new File('*_genomic.fna.gz'),
  output: new File('*_genomic.fna')
}, ({ input }) => shell(`bgzip -d ${input}`) )

const mpileupAndCall = task({
  input: [new File('*_genomic.fna'), new File('*.bam'), new File('*.bam.bai')],
  output: new File('variants.vcf')
}, ({ input }) => shellPipe([
  `samtools mpileup -uf ${input[0]} ${input[1]}`,
  'bcftools call -c - > variants.vcf'
]) )

const pipeline = parallel({
  taskLists: [[samples, fastqDump], [downloadReference, bwaIndex]],
  next: join(alignAndSort, samtoolsIndex, decompressReference, mpileupAndCall)
})
pipeline()
```

### DSL Support

There is early prototypic DSL support using YAML syntax:

```yaml
fastqDump:
  input: "**/*.sra"
  output:
    - "*_1.fastq.gz"
    - "*_2.fastq.gz"
  task: "`fastq-dump --split-files --skip-technical --gzip ${input}`"
```

And consumed with:

```js
const defs = yaml.safeLoad(fs.readFileSync('pipeline.yaml', 'utf-8'))

const fastqDump = taskYAML(defs.fastqDump)
```

### API

**task**: first param is an object with `input` and `ouput`, second is a callack
function that will be passed a resolved version of the first param object

**join**(...tasks): given an arbitrary number of tasks, return a new task

**parallel**([lists], next): given an array of arrays of tasks, and a next task
run the lists tasks in parallel then call next. returns a task.
