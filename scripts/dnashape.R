#!/usr/bin/env Rscript
# Analysis F2: DNAshapeR on one FASTA; writes <fasta>.MGW/.Roll/.HelT/.ProT next to it.
# Called by src/centromere/grammar.py::run_dnashape.
args <- commandArgs(trailingOnly = TRUE)
suppressPackageStartupMessages(library(DNAshapeR))
invisible(getShape(args[1], shapeType = c("MGW", "Roll", "HelT", "ProT")))
