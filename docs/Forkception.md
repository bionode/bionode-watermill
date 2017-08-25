# Forkception

Fork is a special case within bionode-watermill because it requires to know 
everything that proceeds the tasks within the `fork`. In current API there is
 no easy way to get graph shape before actually running the `pipeline.js` 
 scripts. Therefore, `fork` uses `join` to construct independent lineages for
  each task within the `fork`.