const ncbi = require('bionode-ncbi')
const wrapper = require('bionode-wrapper')
const waterwheel = require('bionode-waterwheel')

const { task, join, run } = waterwheel
const { stdout, stdin, file, directory } = waterwheel.types

// Can be passed in from CLI
// params for final pipeline call. params are things that do not
// change how items are passed between processes, but decide output for the
// pipeline as a whole. For example, species name and reads accession.
const pipelineParams = {
  specie: 'Salmonella-enterica',
  readsID: '2492428',
  output: 'Salmonella-enterica.vcf'
}

const sra = task({
  output: stdout()
}, ncbi.download('sra', '{params.readsID}'))

const reference = task({
  output: stdout()
}, ncbi.download('assembly', '{params.specie}'))
const bwa_index = task({
  input: file(),
  output: file()
}, wrapper('bwa index {input}'))


// some tools you cannot decide output file
// so, tell waterwheel to stream a file as stdout
// there is --stdout for fastq-dump to give a streamed fastq, but trimmomatic
// wants reads_1 and reads_2
// tools that are not bionode will need to wrapper()ed
const extract = task({
  input: file('reads.sra'),
  output: file([1, 2].map(n => `reads_${n}.fastq.gz`))
}, wrapper('fastq-dump --split-files --skip-technical --gzip {input}'))

const trim = task({
  input: file([1, 2].map(n => `reads_${n}.fastq.gz`)),
  output: {
    pe: file([1, 2].map(n => `reads_${n}.trim.pe.fastq.gz`),
    se: file([1, 2].map(n => `reads_${n}.trim.se.fastq.gz`)
  },
  opts: {
    adapters: '../adapters'
  }
}, wrapper('''
  trimmomatic PE -phred33 \
  {input[0]} {input[1]} \
  {output.pe[0]} {output.se[0]} \
  {output.pe[1]} {output.se[1]} \
  ILLUMINACLIP:{opts.adapters}/TruSeq3-PE.fa:2:30:10 \
  LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36 \
'''))

const merge = task({
  input: file(),
  output: stdout()
}, wrapper('seqtk mergepe {input}'))

const gzip = task({
  input: stdin(),
  output: stdout()
}, wrapper('gzip - > {output}'))

// Branching: these two filtering types produce the same type of output
// So we can pass an array, and define them both under the same input/output
const filter = task({
  input: file()
  output: file()
  opts: {
    tmpDir: directory()
  }
}, [wrapper('''
  kmc -k{params.KMERSIZE} -m{params.MEMORYGB} -t{params.THREADS} {input} reads.trim.pe.kmc {opts.tmp} \
  kmc_tools filter reads.trim.pe.kmc -cx{params.MINCOVERAGE} {input} -ci0 -cx0 {output} \
'''), wrapper('''
  load-into-counting.py -N 4 -k {params.KMERSIZE} -M {params.MEMORYGB}e9 -T {params.THREADS} reads.trim.pe.fastq.gz.kh {input} \
  abundance-dist.py reads.trim.pe.fastq.gz.kh {input} reads.trim.pe.fastq.gz.kh.hist \
  filter-abund.py -T {params.THREADS} -C ${MINCOVERAGE} reads.trim.pe.fastq.gz.kh -o {output} {input} \
''')])

// file() dependencies will stop stream.
// for example, need to wait on an index file to be made before aligning
const bwa_mem = task({
  input: {
    reference: file(),
    index: file(), // check for reference indexing
    sample: file() // output of filter
  }
  output: stdout()
}, wrapper('bwa mem {input.reference} {input.sample}'))

// This is where streams in Node can really show
// In snakemake or Nextflow, this would be bwa mem | samtools view
// Which is less modular, reproducible, containerizable
const samtools_view = task({
  input: stdin(),
  output: stdout()
}, wrapper('samtools view -Sbh {input}'))

const samtools_sort = task({
  input: stdin(),
  output: stdout()
}, wrapper('samtools sort {input}'))

const samtools_index = task({
  input: stdin(),
  output: stdout()
}, wrapper('samtools index {input}'))

// these will all get piped through each other:
// bwa_mem().pipe(samtools_view()).pipe(samtools_sort()).pipe(samtools_index())
const align = join([bwa_mem, samtools_view, samtools_sort, samtools_index])

const samtools_mpileup = task({
  input: {
    bam: file(),
    index: file(),
    reference: file()
  },
  output: stdout()
}, wrapper('samtools mpileup -uf {input.bam} {input.reference}'))

const bcftools_call = task({
  input: stdin(),
  output: stdout()
}, wrapper('bcftools call -c {input}'))

const callVariants = join([samtools_mpileup, bcftools_call])

// need a way to use output of another task as input for this one
const pipeline = join(
  [sra, extract, trim, merge],
  // this creates a branching
  filter,
  align,
  callVariants
)

// Run the whole pipeline, passing in params
run(pipelineParams, pipeline).pipe(task(fs.createWriteStream('{output}')))
