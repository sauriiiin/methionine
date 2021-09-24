

##### CARBON SOURCE
data.cbn <- data
data.sum.cbn <- data.sum
strain.labs.cbn <- strain.labs2

save(data.cbn, data.sum.cbn, strain.labs.cbn,
     file = 'figures/final/data.RData')

##### JAKES MUTANTS
data.jm <- data
data.sum.jm <- data.es
strain.labs.jm <- strain.labs2

load(file = 'figures/final/data.RData')
save(data.cbn, data.sum.cbn, strain.labs.cbn,
     data.jm, data.sum.jm, strain.labs.jm,
     file = 'figures/final/data.RData')

##### BISMUTH BiGGY
strain.labs.bis <- strain.labs2

load(file = 'figures/final/data.RData')
save(data.cbn, data.sum.cbn, strain.labs.cbn,
     data.jm, data.sum.jm, strain.labs.jm,
     data.bis, strain.labs.bis,
     file = 'figures/final/data.RData')

##### RESPIRATION
data.res <- data
data.res.gc <- gc.pred 
data.res.gc.sum <- gc.res.sum[,-c(11:16)]
data.res.gc.sum2 <- gc.res.sum2
strain.labs.res <- strain.labs2

load(file = 'figures/final/data.RData')
save(data.cbn, data.sum.cbn, strain.labs.cbn,
     data.jm, data.sum.jm, strain.labs.jm,
     data.bis, strain.labs.bis,
     data.res, data.res.gc, data.res.gc.sum, data.res.gc.sum2, strain.labs.res,
     file = 'figures/final/data.RData')

##### HOMOCYSTEIN - GINA
data.bioc <- data.gina
data.bioc.ttest <- ttest.res

load(file = 'figures/final/data.RData')
save(data.cbn, data.sum.cbn, strain.labs.cbn,
     data.jm, data.sum.jm, strain.labs.jm,
     data.bis, strain.labs.bis,
     data.res, data.res.gc, data.res.gc.sum, data.res.gc.sum2, strain.labs.res,
     data.bioc, data.bioc.ttest,
     file = 'figures/final/data.RData')

##### DELETION SCREEN
data.del.diff <- data.diff
data.del.tc <- data.som2
data.del.diff.dist <- diff.dist2
strain.labs.del <- strain.lab.del

load(file = 'figures/final/data.RData')
save(data.cbn, data.sum.cbn, strain.labs.cbn,
     data.jm, data.sum.jm, strain.labs.jm,
     data.bis, strain.labs.bis,
     data.res, data.res.gc, data.res.gc.sum, data.res.gc.sum2, strain.labs.res,
     data.bioc, data.bioc.ttest,
     data.del.diff, data.del.diff.dist, data.del.tc, strain.labs.del,
     file = 'figures/final/data.RData')

##### PLASMID VALIDATION
data.pv <- data
data.sum.pv <- data.sum

load(file = 'figures/final/data.RData')
save(data.cbn, data.sum.cbn, strain.labs.cbn,
     data.jm, data.sum.jm, strain.labs.jm,
     data.bis, strain.labs.bis,
     data.res, data.res.gc, data.res.gc.sum, data.res.gc.sum2, strain.labs.res,
     data.bioc, data.bioc.ttest,
     data.del.diff, data.del.diff.dist, data.del.tc, strain.labs.del,
     data.pv, data.sum.pv,
     file = 'figures/final/data.RData')

##### MODELING
load(file = 'figures/final/data.RData')
save(data.cbn, data.sum.cbn, strain.labs.cbn,
     data.jm, data.sum.jm, strain.labs.jm,
     data.bis, strain.labs.bis,
     data.res, data.res.gc, data.res.gc.sum, data.res.gc.sum2, strain.labs.res,
     data.bioc, data.bioc.ttest,
     data.del.diff, data.del.diff.dist, data.del.tc, strain.labs.del,
     data.pv, data.sum.pv,
     data.mdl,
     file = 'figures/final/data.RData')

