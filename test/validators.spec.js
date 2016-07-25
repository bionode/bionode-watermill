'use strict'

const { assert } = require('chai')
const path = require('path')
const intoStream = require('into-stream')
const fs = require('fs')

const task = require('../lib/Task.js')
const validators = require('../lib/validators.js')

const existsAndNonEmptyFile = path.resolve(__dirname, 'exists.txt')

describe('validators', function() {
  describe('existence check', function() {
    it('should resolve if the specific file exists', function() {
      let myFile = { file: existsAndNonEmptyFile }
      const result = validators.existenceCheck(myFile)

      assert.isOk(myFile.resolved, 'file did not resolve as expected')
      assert.isOk(result)
    })

    it('should not resolve if the file does not exist', function() {
      const myFile = { file: 'nowhere-to-be-found.txt' }
      const result = validators.existenceCheck(myFile)

      assert.isNotOk(myFile.resolved, 'non-existent file was resolved incorrectly')
      assert.isNotOk(result)
    })

    it('should resolve files matching a pattern', function() {
      let myFiles = { file: '*_foo.txt', globOpts: { cwd: __dirname } }

      const result = validators.existenceCheck(myFiles)

      const shoulds = {
        '1_foo.txt': false,
        '2_foo.txt': false
      }

      myFiles.file.forEach((file) => shoulds[file] = true)

      assert.isOk(Object.keys(shoulds).map(key => shoulds[key]).every(item => item === true), 'Did not match expected files')
      assert.isOk(result)
    })
  })

  describe('null check helper', function() {
    it('should gracefully handle a non-existent file', function() {
      const result = validators.isFileNull('the-future.txt')

      assert.isNotOk(result)
    })

    it('should return false for an empty file', function() {
      const result = validators.isFileNull(path.resolve(__dirname, 'empty.txt'))

      assert.isNotOk(result)
    })

    it('should return true for a non-empty file', function() {
      const result = validators.isFileNull(existsAndNonEmptyFile)

      assert.isOk(result)
    })
  })

  describe('null check', function() {
    it.skip('should gracefully handle a non-resolved file', function() {

    })

    it.skip('should return false for an empty file', function() {

    })

    it.skip('should return true for a non-empty file', function() {

    })
  })
})

const onlyNumbers = ({ file }) => {
  const contents = fs.readFileSync(file[0], 'utf-8')
  for (let line of contents.split('\n')) {
    const code = line.charCodeAt(0)

    // Not 0-9
    if (code < 48 || code > 57 ) {
      return false
    }
  }

  return true
}

describe('custom validators', function() {
  describe('custom file validators', function() {
    it.skip('should invalidate a task', function(done) {
      const task1 = task({
        input: { value: intoStream(['1', '2', 'a', '3'].map(i => i + '\n')) },
        output: { file: 'validators.spec.1.txt' },
        skippable: false,
        validators: {
          file: [ onlyNumbers ]
        }
      }, ({ input }) => input.pipe(fs.createWriteStream('validators.spec.1.txt')) )

      task1().on('finish', function() {
        const output = this.output()

        assert.isNotOk(output.resolved)

        done()
      })
    })

    it.skip('should validate a task', function(done) {
      const task1 = task({
        input: { value: intoStream(['1', '2', '5', '3'].map(i => i + '\n')) },
          output: { file: 'validators.spec.2.txt' },
        skippable: false,
        validators: {
          file: [ onlyNumbers ]
        }
      }, ({ input }) => input.pipe(fs.createWriteStream('validators.spec.2.txt')) )

      task1().on('finish', function() {
        const output = this.output()

        assert.isOk(output.resolved)

        done()
      })
    })

  })
})
