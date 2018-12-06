# make plotly pileup ############################################################################
## file: file path of bigWig file
## range: integer or numberic vector of start and end
##


## annoFile: file path to BED file
## chr: character of chromosome
## range: integer or numeric vector of start and end

makeAnnotationTrack <- function(annoFile, chr, d.range, y.index, col, yaxis) { 
    
    importRange <- GRanges(seqnames = chr, ranges = IRanges(start = d.range[1], end = d.range[2]))
    anno <- readRDS(annoFile)
    anno.gr <- makeGRangesFromDataFrame(anno, keep.extra.columns = TRUE)
    # #names(anno.gr) <- anno.gr$name
    anno.ovl <- as.matrix(findOverlaps(anno.gr, importRange))
    anno.df <- as.data.frame(anno.gr[anno.ovl[,1]])
    anno.df <- anno.df[,c("pk_id", "seqnames", "start", "end", "peak.type", "gene.name")]
    
    
    
    # rect <- list(
    #     type = "rect",
    #     opacity = 1,
    #     line = list(color = col, width = 0.5),
    #     xref = "x",
    #     yref = paste0("y", yaxis)
    # )
    # 
    # rects <- list()
    # num <- nrow(anno.df)
    # for (i in 1:num) {
    #     rect[["x0"]] <- anno.df$start[i]
    #     rect[["x1"]] <- anno.df$end[i]
    #     rect[c("y0")] <- y.index - 0.45
    #     rect[c("y1")] <- y.index + 0.45
    #     
    #     rect[["fillcolor"]] <- col
    #     
    #     rects <- c(rects, list(rect))
    # }
    # 
    # 
    # anno <- list(
    #     xref = "x",
    #     yref = paste0("y", yaxis),
    #     showarrow = FALSE
    # )
    # 
    # annos <- list()
    # for (i in 1:nrow(anno.df)) {
    #     
    #     str <- anno.df$name[i]
    #     
    #     anno[["xanchor"]] <- "middle"
    #     anno[["text"]] <- str
    #     anno[["x"]] <- anno.df$start[i] + anno.df$width[i]
    #     anno[["y"]] <- y.index
    #     annos <- c(annos, list(anno))
    # 
    # }
    
    if(nrow(anno.df) > 0){
        m <- anno.df[,c("start", "end")]
        m$filler <- NA
        g.lines <- as.vector(t(m))
        df <- data.frame(xcoord = g.lines, 
                         ycoord = rep(y.index, length(g.lines)), 
                         text = paste0(rep(anno.df$gene.name, each = 3), "\n", rep(anno.df$pk_id, each = 3))
                        )
    } else {
        df <- data.frame(xcoord = seq(d.range[1], d.range[2]), ycoord = NA)
    }
    
    return(df)
    #return(list( rect = rects, anno = annos))
}
    
