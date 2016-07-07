'use strict'

const { assert } = require('chai')
const path = require('path')

const task = require('../lib/Task.js')
const validators = require('../lib/validators.js')

describe('validators', function() {
  describe('existence check', function() {
    it('should resolve if the specific file exists', function() {
      let myFile = { file: path.resolve(__dirname, 'exists.txt') }
      myFile = validators.existenceCheck(myFile)

      assert.isOk(myFile.resolved, 'file did not resolve as expected')
    })

    it('should not resolve if the file does not exist', function() {
      let myFile = validators.existenceCheck({ file: 'nowhere-to-be-found.txt' })

      assert.isNotOk(myFile.resolved, 'non-existent file was resolved incorrectly')
    })

    it('should resolve files matching a pattern', function() {
      let myFiles = { file: '*_foo.txt', globOpts: { cwd: __dirname } }

      myFiles = validators.existenceCheck(myFiles)

      const shoulds = {
        '1_foo.txt': false,
        '2_foo.txt': false
      }

      myFiles.file.forEach((file) => {
        shoulds[file] = true
      })


      assert.isOk(Object.keys(shoulds).map(key => shoulds[key]).every(item => item === true), 'Did not match expected files')
    })
  }) 
})
