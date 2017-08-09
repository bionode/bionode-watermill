'use strict'

const orchestratorUidGenerator = (tasks) => {
  const multiHash = []
  for (const task in tasks) {
    multiHash.push(tasks[task].info.uid)
  }
  return multiHash
}

module.exports = orchestratorUidGenerator
