'use strict'

const chai = require('chai')
const chaiAsPromised = require('chai-as-promised')
chai.use(chaiAsPromised)
const { assert, should } = chai
should()

const applicator = require('../lib/utils/applicator.js')

const tuc = (str) => new Promise(resolve => resolve(str.toUpperCase()))

describe('Applicator', function() {
  it('should resolve a single string', () =>
    applicator('foo', tuc).should.eventually.equal('FOO')
  )

  it('should resolve an array of strings', () =>
    applicator([ 'a', 'b' ], tuc).should.eventually.deep.equal([ 'A', 'B' ])
  )

  it('should handle plain objects', () =>
    applicator({ a: 'a', b: 'b' }, tuc).should.eventually.deep.equal({ a: 'A', b: 'B' })
  )

  it('should handle arrays made of objects', () =>
     applicator([
       {a: 'a'},
       {b: 'b'}
     ], tuc).should.eventually.deep.equal([
       {a: 'A'},
       {b: 'B'}
     ])
  )

  it('should reject if one of the items fails', () =>
    applicator(['a', 1], tuc).should.be.rejected
  )

  it('should handle objects with values as arrays', () =>
     applicator({
       letters: ['a', 'b', 'c']
     }, tuc).should.eventually.deep.equal({
       letters: ['A', 'B', 'C']
     })
  )

  it('should handle a complex mix of nested objects and arrays', () =>
     applicator({
       obj: {
         arrOfObjs: [{
           arr: [ 'a', 'b' ]
         }, {
           str: 'foo'
         }]
       }
     }, tuc).should.eventually.deep.equal({
       obj: {
         arrOfObjs: [{
           arr: [ 'A', 'B' ]
         }, {
           str: 'FOO'
         }]
       }
     })
  )
})
