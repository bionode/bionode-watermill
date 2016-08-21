'use strict'

const crypto = require('crypto')

const hash = (data) => crypto.createHash('sha256').update(data).digest('hex')

module.exports = hash
