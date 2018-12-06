# make plotly pileup ############################################################################
## file: file path of bigWig file
## range: integer or numberic vector of start and end
##

getBwScores <- function (file, d.range, maxbp = 10000, chr) {
    
    importRange <- GRanges(seqnames = chr, ranges = IRanges(start = d.range[1], end = d.range[2]))
    gr.bw <- import(con = file, format = "bigWig", which = importRange)
    
    # if no signal from gr.bw, e.g. goes past chromosome length, then just return zeros
    if (length(gr.bw) > 0){
        
        if ((d.range[2] + 1 - d.range[1])  > maxbp){
            
            binned <- unlist(tile(importRange, n = maxbp))
            ovl <- as.data.frame(findOverlaps(binned, gr.bw))
            ovl$score <- gr.bw$score[ovl$subjectHits]
            ovl.tibble <- ovl %>% group_by(queryHits) %>% summarise( avgScore = mean(score))
            binned$score <- 0
            binned$score[ovl.tibble$queryHits] <- ovl.tibble$avgScore
            x <- mid(ranges(binned))
            y <- binned$score
            
        } else {
            
            cov <- GenomicRanges::coverage(gr.bw, weight = "score")[[chr]]
            x <- seq(d.range[1], d.range[2])
            y <- as.numeric(cov[ranges(importRange)])
        }
        
    } else {
        x <- seq(d.range[1], d.range[2])
        y <- rep( 0, length(x))
    }

    stopifnot(length(y) == length(x))
    print(paste0("length of x and y: ", length(x)))
    return(list(x = x, y = y))
}

makePileupTrack <- function(file, d.range, p.range, col, chr, type) {    
    
    
    
    if (type == "line"){
        xy <- getBwScores(file = file, d.range = d.range, chr = chr)
        x.plot <- xy$x
        y.plot <- xy$y
        y.mx <- round(max(y.plot[which(x.plot > p.range[1] & x.plot < p.range[2])]))
        print(paste0("plot, ymax is ", y.mx))
        p <- plot_ly(x = x.plot, y = y.plot, type = "scatter", mode = 'lines',
                     fill = "tozeroy", alpha = 1, hoverinfo = "x", name = paste0("plot ", 1),
                     fillcolor = col, line = list(color = col)
        ) %>% toWebGL()
    } else if (type == "bar"){
        xy <- getBwScores(file = file, d.range = d.range, chr = chr, maxbp = 10000)
        x.plot <- xy$x
        y.plot <- xy$y
        y.mx <- round(max(y.plot[which(x.plot > p.range[1] & x.plot < p.range[2])]))
        print(paste0("plot, ymax is ", y.mx))
        p <- plot_ly(x = x.plot, y = y.plot, type = "bar", alpha = 1,
                     hoverinfo = "x", name = paste0("plot ", 1),
                     marker = list(color = col)
        ) %>% layout(bargap = 0)
    }
    
    return(list(plot=p, ymax=y.mx))

}