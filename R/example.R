pompExample <- function (example, ..., show = FALSE, envir = .GlobalEnv) {
    example <- as.character(substitute(example))
    ep <- paste0("in ",sQuote("pompExample"),": ")
    pomp.dir <- system.file("examples",package="pomp")
    exampleDirs <- getOption("pomp.examples",default=pomp.dir)
    names(exampleDirs) <- exampleDirs
    show <- as.logical(show)
    if (example=="") {
        avlbl <- lapply(exampleDirs,list.files,pattern=".+?R$")
        avlbl <- lapply(avlbl,function(x) gsub("\\.R$","",x))
        for (dir in exampleDirs) {
            cat("examples in ",dir,":\n",sep="")
            print(avlbl[[dir]])
        }
    } else {
        evalEnv <- list2env(list(...))
        file <- c(lapply(exampleDirs,list.files,
                         pattern=paste0(example,".R"),
                         full.names=TRUE),
                  recursive=TRUE)
        if (length(file)<1) {
            stop(ep,"cannot find file ",
                 sQuote(paste0(example,".R")),call.=FALSE)
        }
        if (length(file)>1) {
            warning(ep,"using ",sQuote(file[1])," from ",sQuote(names(file)[1]),call.=FALSE)
        }
        if (show) {
            file.show(file[1])
            return(invisible(NULL))
        }
        objs <- source(file[1],local=evalEnv)
        if (is.null(envir)) {
            obj <- setNames(lapply(objs$value,get,envir=evalEnv),objs$value)
        } else if (is.environment(envir)) {
            for (i in seq_along(objs$value)) {
                assign(objs$value[i],
                       get(objs$value[i],envir=evalEnv),
                       envir=envir)
            }
            cat("newly created object(s):\n",objs$value,"\n")
            obj <- NULL
        } else {
            stop(ep,sQuote("envir")," must be an environment or NULL",call.=FALSE)
        }
        invisible(obj)
    }
}
