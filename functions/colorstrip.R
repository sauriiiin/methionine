
colorstrip <- function(p, c) {
  g <- ggplot_gtable(ggplot_build(p))
  stripr <- which(grepl('strip-t', g$layout$name))
  fills <- c
  k <- 1
  for (i in stripr) {
    if (length(which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))) != 0) {
      j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
      g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
      k <- k+1
    }
  }
 return(g)
}
