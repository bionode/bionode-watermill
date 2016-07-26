'use strict'

const Promise = require('bluebird')
const path = require('path')

const chai = require('chai')
const chaiAsPromised = require('chai-as-promised')
chai.use(chaiAsPromised)
const { assert, should } = chai
should()

const twd = path.resolve(__dirname, 'files', 'is-non-null')

const isNonNull = require('../lib/validators/is-non-null.js')

describe('is-non-null', function() {
  // Passes locally but fails on Travis
  it('Should resolve if file has contents', () =>
    isNonNull(path.resolve(twd, 'non-empty.txt')).should.be.fulfilled
  )

  it('Should reject if file is null', () =>
    isNonNull(path.resolve(twd, 'empty.txt')).should.be.rejected
  )

  it('Should reject if the file does not exist', () =>
    isNonNull('nonexistent.txt').should.be.rejected
  )
})
