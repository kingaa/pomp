examples <- c("euler.sir","ou2","rw2","verhulst","dacca")

for (ex in examples) {
  source(file=paste(ex,"R",sep="."))
  save(list=ex,file=paste(ex,"rda",sep="."))
}
