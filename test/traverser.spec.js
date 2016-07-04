'use strict'

const { assert, expect } = require('chai')
const traverser = require('../lib/traverser.js')


const onlyResolveValue = (key, obj) => {
  if (key === 'value') {
    obj = obj[key]
    return obj
  } else {
    return false
  }
}

describe('traverser', function() {
  it('should handle base case for value', function() {
    const obj = traverser({ value: 'foo' }, onlyResolveValue)

    assert(obj === 'foo')
  })

  it('should handle an array', function() {
    const obj = traverser([{
      value: 'foo'
    }, {
      value: 'bar'
    }], onlyResolveValue)

    assert(obj[0] === 'foo')
    assert(obj[1] === 'bar')
  })

  it('should handle an object', function() {
    const obj = traverser({
      one: { value: 'foo' }
    }, onlyResolveValue)

    assert(obj.one === 'foo')
  })

  it('should handle a complicated mix', function() {
    const obj = traverser({
      a: { value: 'a' },
      b: [{
        value: 'b1'
      }, {
        c: [{
          value: 'c1'
        }]
      }]
    }, onlyResolveValue)

    assert.deepEqual(obj, {
      a: 'a',
      b: [
        'b1',
        {
          c: ['c1']
        }
      ]
    })

  })
})
