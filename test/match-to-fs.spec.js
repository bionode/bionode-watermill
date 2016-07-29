'use strict'

const Promise = require('bluebird')
const fs = Promise.promisifyAll(require('fs'))
const path = require('path')

const chai = require('chai')
const chaiAsPromised = require('chai-as-promised')
chai.use(chaiAsPromised)
const { assert, should } = chai
should()

const twd = path.resolve(__dirname, 'files', 'match-to-fs')

const matchToFs = require('../lib/matchers/match-to-fs.js')
const applicator = require('../lib/utils/applicator.js')

describe('match-to-fs', function() {
  it('should resolve a single pattern', () =>
    applicator(
      'foo*',
      (item) => matchToFs(item, 1, { cwd: twd })
    ).should.eventually.equal(
      twd + '/' + 'foo.txt'
    )
  )

  it('should reject if there is no match', () =>
    applicator(
      'nonexistent.txt',
      (item) => matchToFs(item, 1, { cwd: twd })
    ).should.be.rejected
  )

  it('should reject if there are an unexpected amount of matches', () =>
     applicator(
       '*_two.txt',
       (item) => matchToFs(item, 1, { cwd: twd })
     ).should.be.rejected
  )
})
