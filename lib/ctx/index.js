'use strict'

const _ = require('lodash')

exports.mergeCtx = (type) => (currentCtx, results) => {
  const { uid, params } = results
  let newTrajection = []

  // If it was a task (not a join or parallel)
  if (!_.isUndefined(uid) && !_.isUndefined(params)) {
    const paramsString = JSON.stringify(params)
    if (paramsString !== '{}') {
      newTrajection = [paramsString, uid]
    } else {
      newTrajection = [uid]
    }
  }

  let newTrajectory = results.context.trajectory.concat(newTrajection)
  if (type === 'junction') {
    newTrajectory = currentCtx.trajectory.concat(newTrajectory)
  }

  const newContext = {
    type,
    trajectory: newTrajectory
  }

  return newContext
}
