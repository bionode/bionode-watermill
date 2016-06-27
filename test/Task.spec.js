const { assert } = require('chai')
const isStream = require('isstream')

const Task = require('../lib/Task.js')

describe('Task', function() {
  // describe('checking parameters', function(done) {
  //   it('should error if not provided with a function', function(done) {
  //     done(new Error('TODO'))
  //   })
  // })

  describe('provided with a Promise', function() {
    it('should return a stream', function(done) {
      const task = Task(
        {
          input: 100
        },
        ({ input }) => new Promise((resolve, reject) => {
          setTimeout(() => resolve('DATUMS'), input)
        })
      )()

      task
        .on('data', (data) => {
          console.log('data: ' + data)
          done()
        })
        .on('end', (data) => console.log(data))

      assert(isStream(task))

      // done(new Error('TODO'))
    })
  })
})
