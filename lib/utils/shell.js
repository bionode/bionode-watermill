'use strict'

const { spawn } = require('child_process')
const { tab } = require('../utils')

const shell = (cmd, opts = {}, logger) => {
  logger.emit('log', 'Starting: ' + cmd)
  cmd = cmd.split(' ')

  opts = Object.assign({ shell:true }, opts)

  const myProcess = spawn(cmd[0], cmd.slice(1), opts)

  myProcess.on('error', onError)
  myProcess.stdout.on('data', onStdout)
  myProcess.stderr.on('data', onStderr)
  myProcess.on('close', onClose)
  myProcess.on('disconnect', onDisconnect)
  myProcess.on('exit', onExit)

  function onStdout(chunk) {
    for (let line of chunk.toString().split('\n')) {
      logger.emit('log', 'stdout: ' + line, 2)
    }
  }

  function onStderr(chunk) {
    for (let line of chunk.toString().split('\n')) {
      logger.emit('log', 'stderr: ' + line, 2)
    }
  }

  function onClose(code) {
    // console.log('CP closed: ' + code)
    if (code !== 0) {
      logger.emit('log', `Child process exiting with code ${code}`)
      logger.emit('log', 'Exiting...')
      setTimeout(() => process.exit(code), 5000)
    }
  }

  function onExit(code) {
    // console.log('CP exited: ', code)
  }

  function onError(err) {
    logger.emit('log', `Child process error: ${err}`)
  }

  function onDisconnect() {
    logger.emit('log', 'Child process disconnect')
  }

  return myProcess
}

exports.shell = shell

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
