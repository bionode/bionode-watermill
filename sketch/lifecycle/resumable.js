'use strict'

const resumable = (taskState) => taskState.resumbable === 'on' ? true : false

module.exports = resumable
