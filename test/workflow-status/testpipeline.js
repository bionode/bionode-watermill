'use strict'

const { task, join } = require('../..')
const { assert } = require('chai')
const mocha = require('mocha')

const { store } = require('../..')



const anotherTask = task({
  name: 'Run a JS file',
  output: 'ran.txt'
}, () => `node /home/evoxtorm/Desktop/Bionode-watermill/bionode-watermill/test/workflow-status/wait-error.js && touch ran.txt`)

//anotherTask()

// Testing the work-flow Status

describe('Failed pipeline', function () {
  it('should have permanentFailure', (done) =>
    anotherTask()
     .then(() => done(new Error('Should not have succeeded'))
    .catch((err) => {
       try {
          expect(store.getState().workflowState).to.equal('PERMANENT_FAILURE')
          done()
     } catch (err) {
        done(err)
     }
   })
)
).timeout(10000)
})


	// .then(console.log)
	// .catch(console.error)



	
	


	