## I will be testing the ET in BioCro
## Date: Dec 17 2014
library(BioCro)

data(cmi04)
res <- BioGro(cmi04)

cumET <- sum(res$CanopyTrans + res$SoilEvaporation) * 0.1

if(cumET < 200 || cumET > 1000) stop("ET in BioGro is incorrect")

maxStem <- max(res$Stem)
maxLeaf <- max(res$Leaf)
maxRoot <- max(res$Root)
maxRhizome <- max(res$Rhizome)

if(maxStem < 5) stop("maxStem in BioGro is incorrect")
if(maxLeaf < 1) stop("maxLeaf in BioGro is incorrect")
if(maxRoot < 0) stop("maxRoot in BioGro is incorrect")
if(maxRhizome < 2) stop("maxRhizome in BioGro is incorrect")
