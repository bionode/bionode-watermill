'use strict'

const hackfile = require('hackfile')
const fs = require('fs')
const util = require('util')

const { File } = require('./lib/types.js')

const fileName = process.argv[2]

const parsed = hackfile(fs.readFileSync(fileName, 'utf-8'))

console.log(util.inspect(parsed, { showHidden: false, depth: null }))

console.log(eval(parsed[0][1][0][1][0]))
