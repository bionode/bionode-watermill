# Graph visualization

Currently bionode-watermill has a tool that allows the user to visualize the 
pipeline Directed Acyclic Graph (DAG). This tool is available while running 
the scripts in `localhost:8084`. Visualization script is available in 
`viz/index.html`.

## Enabling visualization

Just paste `localhost:8084` on your browser URL.

## Canceling

Now, this tool hangs each script at its end and in order to exit it, just 
hit: `ctrl + c` in the shell terminal where the pipeline script is running. 
Hitting `ctrl + c` will close the connection between the browser and 
`lib/reducers/collections.js`.

## How it works

This basically uses npm package `socket.io` on `lib/reducers/collection.js` 
where graph object is read. This graph object follows graphson-like structure:

```javascript
{
  "graph": {
    "mode": "NORMAL",
    "vertices": [
      //...Array of vertices
    ],
    "edges": [
      //...Array of edges
    ] 
  }
}
```

Vertices represent each task on the current graph visualization (as blue 
circles) and edges represent the flow between tasks (with arrows), i.e., 
which tasks runs after other tasks attending to the orchestration shape.
For instance, if we have `join(task1, task2, task3)`, this will render 3 
vertices (blue circles) and 2 edges (arrows) like this:

```javascript
task1 --> task2 --> task3
```

## Node hovering

On node hovering it is possible to obtain metadata information on tasks 
`props`, such as:

- minUid: The short `uid` that matches the filesystem under `./data/<uid>`
- Kind: Currently this differentiates if it is a `task` or a `junction`, 
because `junction`s are added as a orange circle with no relevant proprieties, 
just to make the required connection on the graph.
- Task name: The name given by the user under `props` object while defining a
 task.
- Output(s): The resolved outputs (`resolvedOutput`) from the task run.
- Input(s): The resolved inputs (`resolvedInput`) from the task run.
- params: The `params` given to task `props`
- Output folder: The output folder relative path.
- Cmd: The actual operation passed to shell (when passed) or the 
`operationCreator` function.


## Limitations

* This visualization tool does not represent input/output flow, just
 the pipeline flow.
 * It is also unable to know what was ran from what is running and from 
 what is failing.
