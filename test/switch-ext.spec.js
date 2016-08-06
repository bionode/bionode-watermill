'use strict'

const { assert } = require('chai')

const switchExt = require('../lib/utils/switch-ext.js')

describe('Switch Extension', function() {
  it('Should switch file extension', function() {
    assert.equal(switchExt('foo.txt', 'bar'), 'foo.bar')
  })

  it('Should throw if bad params are passed')
  it('Should handle the regex not being matched')
  it('should handle multiple dots')
})
