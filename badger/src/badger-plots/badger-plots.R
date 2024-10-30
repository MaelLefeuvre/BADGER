#!/usr/bin/env Rscript

if (sys.nframe() == 0){
  library(badger.plots)
  .results <- badger.plots::main()
}