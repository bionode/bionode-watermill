# Quick Start

If you are one of those people that just want to get a quick idea of 
what bionode-watermill can do for you, you have just come to the right place! 

## Example `task`

```javascript
const myTask = task({
  input: '*.fas',  //fasta files present in your current working directory
  output: '*.fas'
}, ({ input }) => {
  const stringInputs = input.join(' ') //puts all files in input array in a
  // single string
  return `cat ${stringInputs} > concatenated.fas`
}  
)
```

This is a task that concatenates (using Unix `cat`) the input files into a 
single output file.

## Example _pipeline_

```javascript
// some tasks are defined before defining the pipeline
const pipeline = join(
  fetchInputs, // a task that gathers input fastas
  junction(
    myTask, // our previous task
    blast   // a task that performes blast of all fastas
  ),
  reportedFastas,  // a task that makes a nice report taking as inputs the
  // blast results and the concatenated fasta in a single file
  fork(
    splitFastas,  // a task that splits the concatenated fastas into its 
    // original fastas
    splitReportedFastas // a task that splits the reported fastas with blast 
    // results into separated reports
    )
  otherCoolAnalysis // finally a task that does some other cool analysis
)
```

Diagram of the code above:

![](https://github.com/bionode/bionode-watermill/blob/docs_update/docs/imgs/quick-start-img.png)

In the above represented `pipeline`, first bionode-watermill will fetch input
 files using some method that you see fit, then the pipeline will split into
  two
  branches. The first branch will perform `myTask`, while the second branch 
  will perform task `blast` and then the pipeline gets merged again into a 
  single
   branch in `reportedFastas` to make that nice report of yours (**that is what 
   `junction` does!**). 
   
   Then the pipeline will split again into two branches and will 
   perform two other tasks `splitFastas` and `splitReportedFastas` but this 
   time it will not get merged into a single branch again. Instead it will 
   duplicate the following task `otherCoolAnalysis` for each one of the 
   branches (**that what `fork` does!**). 
   
   But wait! What assures that this code
    is really run in the proper order? **That is `join`!** Notice how `join` 
    is wrapping everything else! It is in fact making sure everything is run 
    in the order you desire. **Commas are like pauses in join**, it makes 
    sure the pipeline waits for the task before to be finished.

## Running a pipeline

After having making your own script named `pipeline.js` (or whatever you want
 the name to be), just run the following command on your Unix terminal.

```shell
node pipeline.js
```

> There some other logging options that you can explore [here](Loggers.md).