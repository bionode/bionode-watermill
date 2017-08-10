'use strict'

const orchestratorUidGenerator = (tasks) => {
  return tasks.map(task => task.info.uid)
}

module.exports = orchestratorUidGenerator
