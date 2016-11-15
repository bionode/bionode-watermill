'use strict'

const { task } = require('bionode-watermill')

// Generic trimmomatic task does not care whether it receives single or paired
// end reads. It will delegate to a specific task to handle each case. OR maybe
// a bunch of ternaries will do:
// ${isSingle ? '--single' : '--double'}
// probably depends on how task params are organized, ternary might or might not
// be more messy than two delegated tasks.


// In general, we give trimmomatic some reads (whether single or paired), and
// expect back the trimmed equivalents (identified with a .trim. in the filename)
// The adapters to use depends on single vs. paired

// The paired end case:
// input.reads is an array of reads
// output.reads is an array of trimmed reads
// (could have just used "output", but using "output.reads" maintains same
// schema)

// TODO ways to generalize this, (down the line)
// right now, better to be super explicit and deal with consequences of handling
// .fq vs. .fastq, find patterns from that, then move it into core.
const readsFiletype = 'fastq.gz'
const trimmedPESuffix = 'trim.pe'

const pairedTrimmomatic = task(
  {
    input: {
      reads: [1, 2].map(n => `*_${n}.${readsFiletype}`) // readabilty ruined as opposed to `_${n}.fastq.gz`
    },
    // way 1: not using same shape as input. use cases?? better/worse??
    output: [1, 2].map(n => `*_${n}.${trimmedPESuffix}.${readsFiletype}`),
    params: {
      // "hardcode" in params specific to paired trimmomatic here.
      // could use a higher order task (i.e. function that returns a task) to do
      // this also.
      trimmomaticExecutable: path.resolve(__dirname, 'bin', 'trimmomatic-0.36.jar'),
      adapters: path.resolve(__dirname, 'adapters'),
      // or maybe, since the params points to adapters directory, it then depends
      // on bash declaration to choose files specific to paired vs end.
      // probably better to put those here though.
    },
    name: 'paired-end-trimmomatic' // something for the logs. "task instances" of this "name" will be created and added to DAG
  },
  // now declare function that takes resolved "props" which is input which has
  // patterns mapped to absolute file paths, and params. input keeps same shape
  // just replaces patterns with paths.
  ({ input, params }) => {
    // can do stuff, like a "render" function from React, sometimes we use
    // "props" to decide what to render. i.e. ${props.isLoggedIn ? <Dashboard> : <LoginScreen>}
    // deciding where to put logic for single vs paired is tricky,
    // do it here, or have two tasks that get delegated to by something else?

    // in our case, just written a string of bash
    // will be ran by require('child_process').spawn(bin, args, { shell: true })
    // shell: true just runs it as bash -c "what --you --return --here"
    // by default, watermill uses "shell: true" inside
    // TODO document that fact somewhere.

    // here we "hardcode" output filenames to be something super generic,
    // like: reads_2.trim.pe.fastq.gz
    // this has some binding to the variable trimmedPESuffix above
    // utitlies to handle dynamic suffixes etc would be nice

    // safe to be generic "reads" b/c each task runs isolated in its in
    // directory with inputs brought in via symlinks

    return `
java -jar ${params.trimmomaticExecutable} PE -phred33 \
${input[0]} ${input[1]} \
reads_1.${trimmedPESuffix}.fastq.gz reads_1.trim.se.fastq.gz \
reads_2.trim.pe.fastq.gz reads_2.trim.se.fastq.gz \
ILLUMINACLIP:${params.adapters}/TruSeq3-PE.fa:2:30:10 \
LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36`
  }
)

const singleTrimmomatic = task(
  {
    input: {
      reads: '*_.fastq.gz'
    },
    output: '*_.trim.se.fastq.gz',
    params: {
      trimmomaticExecutable: path.resolve(__dirname, 'bin', 'trimmomatic-0.36.jar'),
      adapters: path.resolve(__dirname, 'adapters'),
    },
    name: 'single-end-trimmomatic' // something for the logs. "task instances" of this "name" will be created and added to DAG
  }, ({ input, params }) => `
java -jar ${params.trimmomaticExecutable} SE -phred33 \
${input.reads} reads.trim.se.fastq.gz \
ILLUMINACLIP:TruSeq3-SE:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36`)

// Flow types would be nice here.
// type Reads = [<String>]
// types for paired or single reads, etc.

const trimmomatic = (mode, reads, props) => {
  if (mode === 'single') {
    return singleTrimmomatic
  } else if (mode === 'paired') {
    return pairedTrimmomatic
  }
}

// so now this module exports a function that returns a single or paired
// task. up to dev to know how to use either.
module.exports = trimmomatic

// Lots of others to go about this, could abstract more things,
// but don't want to end up with leaks -
// `_${num}_${suffix1}_${suffix2}_${filetype}` is generalized but so
// unreadable...

// === OLD VERSION: ASSUMES PAIRED END ===

/**
 * Trimming with trimmomatic
 */
const trim = task({
  input: [1, 2].map(n => `*_${n}.fastq.gz`),
  output: [1, 2].map(n => `*_${n}.trim.pe.fastq.gz`),
  params: {
    trimmomaticExecutable: path.resolve(__dirname, 'bin', 'trimmomatic-0.36.jar'),
    adapters: path.resolve(__dirname, 'adapters')
  },
  name: 'Trimming with trimmomatic'
}, ({ input, params }) => `java -jar ${params.trimmomaticExecutable} PE -phred33 \
${input[0]} ${input[1]} \
reads_1.trim.pe.fastq.gz reads_1.trim.se.fastq.gz \
reads_2.trim.pe.fastq.gz reads_2.trim.se.fastq.gz \
ILLUMINACLIP:${params.adapters}/TruSeq3-PE.fa:2:30:10 \
LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36`)

// Maybe add hooks?
// trim.beforeHooks = [ ... ]
// trim.afterHooks = [ .. ]
// custome task lifecycle execution etc, add custom post validations
// right now, post validation checks for non null of resolved outputs
// could also add unique stream watches. for example, if stdout has "error: x
// not found" but the child process gives a 0 exit code (UGH!), then our custom
// stream watcher can say "no wait, this task actually failed"

// module.exports = trim

// Or maybe export task and some metadata?
// TODO define schema for task module exports??

// module.exports = {
//   task: trim,
//   metadata: {}
// }
