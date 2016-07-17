'use strict'

const { store } = require('.')

const actions = require('./lib/get-actions.js')(store.dispatch)

actions.setConfigValue('resume', 'on')
actions.addTask('foo')
actions.addTask('bar')
