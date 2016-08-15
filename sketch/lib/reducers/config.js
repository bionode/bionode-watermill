'use strict'

const SET_CONFIG_VALUE = 'SET_CONFIG_VALUE'

const defaultState = {
  threads: 1,
  resume: 'on', // force-on, on, off, force-off
  workdir: './data',
  waterfile: './waterfile.json',
  statefile: 'state.json',
  container: null
}

const reducer = (state = defaultState, action) => {
  switch (action.type) {
    case SET_CONFIG_VALUE:
      // TODO handle these being undefined
      const { key, val } = action

      if (!(key in defaultState)) {
        throw new Error('Cannot set a new key in config')
        return state
      }

      const newVal = {}
      newVal[key] = val

      return Object.assign(state, newVal)
      break
    default:
      return state
  }
}

reducer.setConfigValue = (key, val) => ({ type: SET_CONFIG_VALUE, key, val })

// So that task can get default values
reducer.defaultState = defaultState

module.exports = reducer

