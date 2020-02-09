# this is to check the data from Chen, Lu, et al. "Genetic drivers of  
# epigenetic and transcriptional variation in human immune cells." 
# Cell 167.5 (2016): 1398-1414.

data.dir = "~/research/data/EGA"

meta = readRDS(file.path(data.dir, "META2/dat.rds"))
dim(meta)

meta[1:2,1:5]

names(meta)
