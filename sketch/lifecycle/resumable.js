'use strict'

const resumable = (taskState) => new Promise((resolve, reject) => {
  if (!taskState) reject('No taskState was provided to resumable(taskState)')

  const { resume } = taskState

  if (!resume) reject('Task state has no `resume` property')

  resolve(resume === 'on' ? true : false)
})

module.exports = resumable
