"use strict"

const {
  task,
  join
} = require("../..")
const {
  assert
} = require("chai");
const mocha = require("mocha")
const path = require('path')

const {
  store
} = require("../..")
var expect = require("chai").expect

const anotherTask = task({
    name: "Run a JS file",
    output: "ran.txt"
  },
  () =>
  `node ${path.resolve(__dirname, 'wait-error.js')} && touch ran.txt`
)


anotherTask()

//Testing the work-flow Status

describe("Failed pipeline", function() {
  it("Should have permanentFailure", done =>
    anotherTask()
    .then(() => done(new Error("Should not have succeeded")))
    .catch(err => {
      try {
        expect(store.getState().workflowState).to.equal("permanentFailure")
        done();
      } catch (err) {
        done(err)
      }
    })).timeout(6000)
})