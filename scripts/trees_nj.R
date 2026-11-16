#!/usr/bin/env Rscript
# Analysis H step 2.4: NJ trees on the K80 matrices, Robinson-Foulds and Mantel.
# Usage: Rscript scripts/trees_nj.R <trees dir>   (reads dist_cen.tsv, dist_genome.tsv)
# Writes cen_nj.nwk, genome_nj.nwk (midpoint-rooted) and tree_compare.tsv.
suppressPackageStartupMessages({ library(ape); library(vegan) })

args <- commandArgs(trailingOnly = TRUE)
dir <- if (length(args)) args[1] else "output/scrap/trees"

read_dist <- function(f) {
  m <- read.delim(file.path(dir, f), row.names = 1, check.names = FALSE)
  as.matrix(m)
}

# Root at the midpoint of the longest tip-to-tip path (ape has no midpoint rooting).
midpoint <- function(tr) {
  d <- cophenetic(tr)
  ij <- which(d == max(d), arr.ind = TRUE)[1, ]
  half <- max(d) / 2
  path <- nodepath(tr, ij[1], ij[2])
  # walk from tip i towards the root and down to tip j
  acc <- 0
  for (k in seq_len(length(path) - 1)) {
    a <- path[k]; b <- path[k + 1]
    e <- which((tr$edge[, 1] == a & tr$edge[, 2] == b) | (tr$edge[, 1] == b & tr$edge[, 2] == a))
    len <- tr$edge.length[e]
    if (acc + len >= half) break
    acc <- acc + len
  }
  child <- tr$edge[e, 2]                   # node below the edge holding the midpoint
  below <- if (child <= Ntip(tr)) tr$tip.label[child] else extract.clade(tr, child)$tip.label
  # distance from the midpoint to `child`
  to_child <- if (child == a) half - acc else len - (half - acc)
  r <- root(tr, node = if (child <= Ntip(tr)) NULL else child,
            outgroup = if (child <= Ntip(tr)) tr$tip.label[child] else NULL, resolve.root = TRUE)
  root_node <- Ntip(r) + 1
  re <- which(r$edge[, 1] == root_node)
  for (x in re) {
    kid <- r$edge[x, 2]
    tips <- if (kid <= Ntip(r)) r$tip.label[kid] else extract.clade(r, kid)$tip.label
    r$edge.length[x] <- if (setequal(tips, below)) to_child else len - to_child
  }
  r
}

dc <- read_dist("dist_cen.tsv")
dg <- read_dist("dist_genome.tsv")[rownames(dc), colnames(dc)]
tc <- midpoint(nj(as.dist(dc)))
tg <- midpoint(nj(as.dist(dg)))
write.tree(tc, file.path(dir, "cen_nj.nwk"))
write.tree(tg, file.path(dir, "genome_nj.nwk"))

n <- nrow(dc)
rf <- dist.topo(unroot(tc), unroot(tg), method = "PH85")
set.seed(1)
mt <- mantel(as.dist(dc), as.dist(dg), method = "spearman", permutations = 9999)
res <- data.frame(n_taxa = n, rf = as.numeric(rf), rf_max = 2 * (n - 3),
                  nrf = as.numeric(rf) / (2 * (n - 3)), mantel_r = unname(mt$statistic),
                  mantel_p = mt$signif, permutations = mt$permutations)
write.table(res, file.path(dir, "tree_compare.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
