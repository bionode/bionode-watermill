'use strict'

const _ = require('lodash')

exports.mergeCtx = (type) => (currentCtx, results) => {
  const { uid, params } = results

  if (_.isUndefined(uid)) {
    throw new Error('Undefined UID')
  }

  const newTrajection = [uid]

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
