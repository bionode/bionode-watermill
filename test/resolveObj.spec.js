'use strict'

const { assert, expect } = require('chai')
// const resolveObj = require('../lib/resolveObj.js')


describe('resolveObj', function() {
  it.skip('should handle base case for value', function() {
    const obj = resolveObj({ value: 'foo' })

    assert(obj === 'foo')
  })

  it.skip('should handle an array', function() {
    const obj = resolveObj([{
      value: 'foo'
    }, {
      value: 'bar'
    }])

    assert(obj[0] === 'foo')
    assert(obj[1] === 'bar')
  })

  it.skip ('should handle an object', function() {
    const obj = resolveObj({
      one: { value: 'foo' }
    })

    assert(obj.one === 'foo')
  })

  it.skip ('should handle a complicated mix', function() {
    const obj = resolveObj({
      a: { value: 'a' },
      b: [{
        value: 'b1'
      }, {
        c: [{
          value: 'c1'
        }]
      }]
    })

    // console.log(JSON.stringify(obj, null, 2))
  })
})
