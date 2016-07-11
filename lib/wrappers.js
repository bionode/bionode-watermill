'use strict'

const { spawn } = require('child_process')
const { tab } = require('./utils.js')

class Shell {
  constructor(cmd, opts = {}) {
    console.log(tab(1) + 'Starting: ' + cmd)
    cmd = cmd.split(' ')

    this._process = spawn(cmd[0], cmd.slice(1), Object.assign({ shell: true }, opts))

    this._process.on('error', this.onError)
    this._process.stdout.on('data', this.onStdout)
    this._process.stderr.on('data', this.onStderr)
    this._process.on('close', this.onClose)
  }

  onStdout(chunk) {
    for (let line of chunk.toString().split('\n')) {
      console.log(tab(2) + 'stdout: ' + line)
    }
  }

  onStderr(chunk) {
    for (let line of chunk.toString().split('\n')) {
      console.log(tab(2) + 'stderr: ' + line)
    }
  }

  onClose(code) {
    if (code !== 0) {
      console.log(`Child process exited with code ${code}`)
    }
  }

  onError(err) {
    console.log('Child process error: ' + err)
  }

  getSpawned() {
    return this._process
  }
}

const shell = (cmd, opts = {}) => {
  console.log(tab(1) + 'Starting: ' + cmd)
  cmd = cmd.split(' ')

  opts = Object.assign({ shell:true }, opts)

  const process = spawn(cmd[0], cmd.slice(1), opts)

  process.on('error', onError)
  process.stdout.on('data', onStdout)
  process.stderr.on('data', onStderr)
  process.on('close', onClose)

  function onStdout(chunk) {
    for (let line of chunk.toString().split('\n')) {
      console.log(tab(2) + 'stdout: ' + line)
    }
  }

  function onStderr(chunk) {
    for (let line of chunk.toString().split('\n')) {
      console.log(tab(2) + 'stderr: ' + line)
    }
  }

  function onClose(code) {
    // console.log('CP closed: ' + code)
    if (code !== 0) {
      console.log(`Child process exited with code ${code}`)
      console.log('Exiting...')
      process.exit(code)
    }
  }

  function onError(err) {
    console.log('Child process error: ' + err)
  }

  return process
}

exports.shell = shell

// exports.shell = (cmd, opts = {}) => {
//   console.log('Starting: ' + cmd)
//   cmd = cmd.split(' ')
//
//   const process = spawn(cmd[0], cmd.slice(1), Object.assign(opts, { shell: true }))
//   process.on('error', msg => console.log(msg))
//
//   process.stdout.on('data', (data) => console.log(`stdout: ${data}`) )
//
//   process.stderr.on('data', (data) => console.log(`stderr: ${data}`) )
//
//   process.on('close', (code) => console.log(`child process exited with code ${code}`) )
//
//   return process
// }

exports.shellPipe = (cmds) => {
  const processes = cmds.map((cmd, i) => {
    console.log('Starting: ' + cmd)
    cmd = cmd.split(' ')
    // needs shell:true for bcftools call > variant.vcf to work
    const process = spawn(cmd[0], cmd.slice(1), { shell: true })

    process.stderr.on('data', (data) => console.log(`stderr-${i}: ${data}`) )

    process.on('close', (code) => console.log(`child process-${i} exited with code ${code}`) )

    return process
  })

  for (let i=0; i < processes.length-1; i++) {
    const current = processes[i]
    const next = processes[i+1]

    current.stdout.on('data', data => next.stdin.write(data))
    current.on('close', () => next.stdin.end())
  }

  return processes[processes.length-1]
}
