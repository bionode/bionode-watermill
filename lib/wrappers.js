'use strict'

const { spawn } = require('child_process')

exports.shell = (cmd, opts = {}) => {
  console.log('Starting: ' + cmd)
  cmd = cmd.split(' ')

  const process = spawn(cmd[0], cmd.slice(1), Object.assign(opts, { shell: true }))
  process.on('error', msg => console.log(msg))

  process.stdout.on('data', (data) => console.log(`stdout: ${data}`) )

  process.stderr.on('data', (data) => console.log(`stderr: ${data}`) )

  process.on('close', (code) => console.log(`child process exited with code ${code}`) )

  return process
}

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
