'use strict'

const { task } = require('../../..')

const fs = require('fs')
const path = require('path')
const { spawn } = require('child_process')

const script = ({ dir, program }) => (strings, ...values) => {
  let content = ''
  strings.forEach((str, i) => {
    content += str + (values[i] || '')
  })

  const scriptLocation = path.resolve(dir, 'script.' + (() => {
    switch (program) {
      case 'bash': return 'sh'
      case 'python': return 'py'
    }
  })())

  fs.writeFileSync(scriptLocation, `
#!/usr/bin/env ${program}
${content.replace(/\n\ {2}/g, '\n')}
  `.trim())

  fs.chmodSync(scriptLocation, '755')

  const cp = spawn(scriptLocation)

  cp.on('error', (err) => console.log(`${scriptLocation} error:`, err))
  cp.on('close', () => console.log(`${scriptLocation} closed`))
  cp.on('exit', () => console.log(`${scriptLocation} exited`))
  cp.on('disconnect', () => console.log(`${scriptLocation} disconnected`))

  cp.stdout.on('data', (chunk) => console.log(`${scriptLocation} stdout: ${chunk}`))
  cp.stderr.on('data', (chunk) => console.log(`${scriptLocation} stderr: ${chunk}`))

  return cp
}

const myTask = task({
  name: 'my task',
  input: '*.lowercase',
  output: '*.uppercase'
}, ({ input, dir }) => script({ dir, program: 'bash' })`
  input="${input}"
  output="${input.replace(/\.lowercase$/, '.uppercase')}"

  cat $input | awk '{print toupper($0)}' > $output
`)

myTask()
  .then(() => console.log('done'))
  .catch(console.error)
