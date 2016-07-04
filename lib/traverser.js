'use strict'

const traverser = (obj, cb) => {
  const traverse = (obj) => {
    if (obj instanceof Array) {
      return obj.map(item => traverse(item))
    }

    for (let key in obj) {
      const mutation = cb(key, obj)

      if (!mutation) {
        obj[key] = traverse(obj[key])
      } else {
        obj = mutation
      }

    }
    
    return obj
  }

  return traverse(obj)
}

module.exports = traverser
