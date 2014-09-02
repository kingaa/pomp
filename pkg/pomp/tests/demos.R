library(pomp)

pdf.options(useDingbats=FALSE)
pdf(file="demos.pdf")

set.seed(47575684)

demos <- list.files(path=system.file("demo",package="pomp"),pattern=".\\.R$",full.names=TRUE)

for (d in demos) {
  source(d,local=TRUE,echo=TRUE)
}

dev.off()
