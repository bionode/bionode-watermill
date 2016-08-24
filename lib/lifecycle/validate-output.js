'use strict'

const {
  isNonNull,
  matchHash,
  matchTime
} = require('../validators')
const applicator = require('../utils/applicator.js')

const { CHECKING_RESUMABLE } = require('../constants/task-status-types.js')

const dispatch = ({ content }) => console.log(content)
const TASK_LOG = 'sdsd'
const SET_VALIDATION = 'sdsds'
const uid = 'sdsd'
const tab = () => ''
const mode = 'sdsds'

const validateOutput = (taskState, logger) => new Promise((resolve, reject) => {
  const { uid, resolvedOutput, status } = taskState
  logger.emit('log', `Running validations for ${uid}`)

  // TODO validate on streams
  if (!taskState.output) {
    logger.emit('log', 'Skipping since no file patterns in output', 1)
    return resolve({
      validated: true,
      validations: {}
    })
  }

  const validArr = [
    applicator(resolvedOutput, isNonNull)
  ]
  // TODO
  if (mode === 'before') {
    validArr.push(applicator(resolvedOutput, matchHash))
    validArr.push(applicator(resolvedOutput, matchTime))
  }
  const validations = Promise.all(validArr)

  // TODO dispatch something to store validations output

  validations
  .then(() => {
    if (status === CHECKING_RESUMABLE) {
      dispatch({
        type: TASK_LOG,
        uid,
        content: tab(1) + 'All validations successful'
      })
    }

    if (mode === 'after') {
      const checksummer = require('../post-validators/checksummer.js')
      return applicator(resolvedOutput, checksummer)
        .then(() => console.log('  ' + 'Wrote output(s) time+hash to watermill.json'))
    }
  })
  .catch((err) => {
    console.log('validations error: ', err)
    logger.emit('log', 'Could not resovle output', 1)
  })
  .then(() => resolve( {
    validated: true,
    validations: {}
  }))
})

module.exports = validateOutput
