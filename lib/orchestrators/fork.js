'use strict'

const _ = require('lodash')
//const Promise = require('bluebird')

const { defaultContext } = require('../constants/default-task-state.js')
const hash = require('../utils/hash.js')
const orchestratorUidGenerator = require('./orchestrator-uid-generator')

const fork = () => {
  return function (...tasks) {
    const forkInvocator = (cb = _.noop, ctx = defaultContext) => {}

    forkInvocator.info = {
      tasks: tasks,  // returns an array of tasks
      type: 'fork',  // used to check the type of orchestrator
      uid: hash(orchestratorUidGenerator(tasks).join(''))
    }

    //console.log("forkInvocator")
    //console.log(forkInvocator.info)

    return forkInvocator
  }

}

module.exports = fork
