'use strict'

const Promise = require('bluebird')
const fs = Promise.promisifyAll(require('fs'))
const checksum = require('checksum')

const checksummer = (path) => new Promise((resolve, reject) => {
  // console.log(tab(1) + 'Now going to generate checksums on resolved files')

  // TODO promisify checksum
  checksum.file(path, (err, sum) => {
    if (err) return reject(err)

    fs.statAsync(path)
      .then((stats) => {
        const objDetails = {
          [path]: {
            hash: sum,
            time: stats.mtime.getTime()
          }
        }

        // TODO promisify
        let currentObj = {}
        try {
          currentObj = JSON.parse(fs.readFileSync('waterwheel.json', 'utf-8'))
        } catch (err) {
          // Probably file does not exist...
          console.log('Assuming b/c nonexistence..')
          // console.log(err)
        }

        const writeableObj = Object.assign({}, currentObj, objDetails)

        fs.writeFileSync('waterwheel.json', JSON.stringify(writeableObj, null, 2))

        console.log(`Wrote info for ${path} to waterwheel.json`)
        resolve()
          // .then(() => {
          //   // console.log(`Wrote info for ${path} to waterwheel.json`)
          //   resolve()
          // })
      })
      .catch(reject)
  })
})

module.exports = checksummer
