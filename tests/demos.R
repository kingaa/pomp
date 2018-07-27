library(pomp)

png(filename="demos-%02d.png",res=100)

set.seed(47575684L)

demos <- list.files(path=system.file("demo",package="pomp"),
                    pattern=".\\.R$",
                    full.names=TRUE)

for (d in demos) {
  source(d,local=TRUE,echo=TRUE)
}

dev.off()
