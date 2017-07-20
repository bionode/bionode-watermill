# bionode-watermill - mappers example pipeline

This example assumes that either you have completed the 
[bionode-watermill tutorial](https://github.com/bionode/bionode-watermill-tutorial)
or that you are already comfortable enough to produce your own pipelines 
using bionode-watermill.

## Requirements for this example

First of all, you will need to clone bionode-watermill if you haven't done it
 yet: `git clone https://github.com/bionode/bionode-watermill`

Also, you need to have installed:

1) **Node.js** version 7 or higher.
2) **[bionode-ncbi](https://github.com/bionode/bionode-ncbi)** (however, this 
one is not really necessary assuming that you 
run the examples from within the bionode-watermill folders).
3) **[sra-toolkit](https://www.ncbi.nlm.nih.gov/books/NBK158900/)**
4) **[bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)**
5) **[bwa](http://bio-bwa.sourceforge.net/)**
6) **[samtools](http://samtools.sourceforge.net/)**

However, don't worry you can use [this](https://github.com/bionode/bionode-watermill-tutorial/tree/master/docker-watermill-tutorial) 
docker image instead. You can also check the [docker hub repo](https://hub.docker.com/r/tiagofilipe12/bionode-watermill-tutorial/) 
if you prefer.

## What for?

This example pipeline intends to download a reference genome of 
_Streptococcus pneumoniae_ and the SRA accession 'ERR045788', and then map the
 reads available in this accession against the reference genome.
 Then, some other processes are suggested but you can configure this pipeline
  according to your own needs (which might be a good starting point to 
  understand the mechanics of bionode-watermill).
  
  This example pipeline is available [here](https://github.com/bionode/bionode-watermill/blob/master/examples/pipelines/two-mappers/pipeline.js).

## The Pipeline
### Fetching sequence data

First, we will define a variable `config` that defines the **name** and
 **URL** of the reference genome of _Streptococcus pneumoniae_. Also this 
 `config` also has the accession of the desired reads (`sraAccession`).
  
```javascript
const config = {
  name: 'Streptococcus pneumoniae',
  sraAccession: 'ERR045788',
  referenceURL: 'http://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/007/045/' +
  'GCF_000007045.1_ASM704v1/GCF_000007045.1_ASM704v1_genomic.fna.gz'
}
```

Then we will download both the reference genome as well as the sample reads.
The first one we download using the `config.referenceURL` previously defined:

```javascript
// first lets get the reference genome for our mapping
const getReference = task({
  params: { url: config.referenceURL },
  output: '*_genomic.fna.gz',
  name: 'Download reference genome for ${config.name}'
}, ({ params, dir }) => {
  const { url } = params
  const outfile = url.split('/').pop()

  // essentially curl -O
  return request(url).pipe(fs.createWriteStream(dir + '/' + outfile))
})
```

While for the second one we will use bionode-ncbi to download the reads:

```javascript
//then get samples to work with
const getSamples = task({
  params: {
    db: 'sra',
    accession: config.sraAccession
  },
  output: '**/*.sra',
  dir: process.cwd(), // Set dir to resolve input/output from
  name: 'Download SRA ${config.sraAccession}'
}, ({ params }) => `bionode-ncbi download ${params.db} ${params.accession}`
)
```

Then since the download from SRA returns a `.sra` file, we have to process with
 `fastq-dump` (sra-toolkit) it in order to get the ` .fastq.gz` files 
 contained in 
 this sra accession.
 
 ```javascript
// extract the samples from fastq.gz
const fastqDump = task({
  input: '**/*.sra',
  output: [1, 2].map(n => `*_${n}.fastq.gz`),
  name: 'fastq-dump **/*.sra'
}, ({ input }) => `fastq-dump --split-files --skip-technical --gzip ${input}`
)
```

Now, note that we are getting reference and samples in parallel using 
**junction** orchestrator:

```javascript
// === PIPELINE ===

const pipeline = join(
  junction(
      getReference,
      join(getSamples,fastqDump)
  ),
  gunzipIt,
  fork(
    join(IndexReferenceBwa, bwaMapper),
    join(indexReferenceBowtie2, bowtieMapper)
  )
)
```

This means that only after both reference and samples have finished 
downloading, and in the case of samples being extracted from `.sra` file, the
 remaining tasks will be executed. In this case the tasks within the **fork**
  orchestrator will be performed after the **junction**.

### Mapping

Now, in fact **fork** isn't really necessary in the above example, however, 
the idea is for you to continue on your own after this example. For instance 
you may try to run `samtools`  after the fork, which will render two 
independent results, one for `bowtie` results and another for `bwa` results. 
The pipeline would look something like this:

```javascript
// === PIPELINE ===

const pipeline = join(
  junction(
      getReference,
      join(getSamples,fastqDump)
  ),
  gunzipIt,
  fork(
    join(IndexReferenceBwa, bwaMapper),
    join(indexReferenceBowtie2, bowtieMapper)
  ),
  samtoolsStuff
)
```

But, back to mapping! In this example we set up mapping in two types of tasks:

* Reference indexing (`IndexReferenceBwa` and `indexReferenceBowtie2`).
* Mapping (`bwaMapper` and `bowtieMapper`).

Then, we 'forced' each indexing to run before each mapping approach by using 
**join** orchestrator. You can check the definition of these tasks 
[here](https://github.com/bionode/bionode-watermill/blob/master/examples/pipelines/two-mappers/pipeline.js#L66-L122).

### The challenge! :boom:

Now that you have tasted a simple pipeline that is able to map sequences from
 NCBI SRA against a reference genome using two different mappers, feel free 
 to test some usages with the output files of the mapping.
 One useful example would be to use samtools to obtain the 'mapping depth' of
  your mapping results. Something like: `samtools faidx` --> `samtools view` 
  --> `samtools sort` --> `samtools index` --> `samtools depth`.


If you run into trouble, you can find an example of this pipeline 
[here](https://github.com/bionode/bionode-watermill/tree/master/examples/pipelines/two-mappers/pipeline_lazy.js).