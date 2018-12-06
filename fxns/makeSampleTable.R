
library(stringr)
library(dplyr)
options(stringsAsFactors = FALSE)

makeSampleTable <- function() {
    data <- list.files("./data", recursive = TRUE)
    data <- grep(".rds$", data, invert = TRUE, value = TRUE)
    st <- as.data.frame(str_split(data, "\\/", simplify = TRUE))[,1:2]
    colnames(st) <- c("experiment", "sample")
    st$bwFile <- data
    st <- st %>% group_by(sample) %>% mutate(replicate = 1:length(bwFile), 
                                             annoFile = grep(".rds", 
                                                             list.files(paste0("./data/",experiment)),
                                                             value = TRUE))
    st <- as.data.frame(st)
    st$bwFile <- paste0("./data/", st$bwFile)
    st$annoFile <- paste0("./data/", st$experiment, "/", st$annoFile)
    st
}