'use strict'

const switchExt = (file, newExt) => {
  const newFile = file.replace(/\.[^.]+$/, '.' + newExt)

  return newFile
}


module.exports = switchExt
