'use strict'

const { task } = require('bionode-watermill')

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

module.exports = trim

// Or maybe export task and some metadata?
// TODO define schema for task module exports??

// module.exports = {
//   task: trim,
//   metadata: {}
// }
