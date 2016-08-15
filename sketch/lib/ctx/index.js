'use strict'

const _ = require('lodash')

exports.defaultCtx = {
  trajectory: []
}

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

  newTrajection = results.trajectory.concat(newTrajection)

  return {
    type,
    trajectory: currentCtx.trajectory.concat(newTrajection)
  }
}
