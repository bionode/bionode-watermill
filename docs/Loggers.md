# Loggers

Currently bionode-watermill has three main logging options.

## Chrome inspector

You can run any standard bionode-watermill script using google chrome 
inspector, like this:

```
node --inspect pipeline.js
```

This will log both to the console and to google chrome. It will give you a 
link where you can see the logs. There are also some google chrome plugins 
that allow for the debugger to open automatically.

## Redux logger

Since bionode-watermill uses `redux` to manage `actions` and `states` and not
 allways what is going on within redux is clear, `redux-logger` greatly helps
  seeing what is happening on each task. You have `previous state`, `action` 
  and `next action` which as it says is useful to access what happens on each
   `reducer` and this might be useful for debugging your pipeline.
   
   This logger is available in bionode-watermill through an environmental 
   variable `REDUX_LOGGER=1`.
   
   ```
REDUX_LOGGER=1 node --inspect pipeline.js
``` 

Here it is strongly recommended to use chrome inspector since `redux-logger` 
prints a lot of outputs and they can be interact with.

## D3 logger or D3 visualization

D3 logger is basically a visualization tool that allows the user to check the
 pipeline shape, outputs/inputs and commands being executed using `d3`. This 
 can be accessed by passing another environmental variable `D3_LOGGER=1`.
 
 ```
D3_LOGGER=1 node pipeline.js
```

This will send `graphson` object to `d3` in order to draw the aspect of the 
pipeline and can be accessed in `http://localhost:8084/`. Other metadata of 
each task will be also shown on hovering (for more details see 
[D3Visualization](D3Visualization.md)).