# Get plotting elements for gene annotations ####################################
## range: integer or numeric vector of start and end
## chr: character of chromosome
## collapseTx: logical. determines which exon reference to retrieve: transcript-level or gene-level

# TO FIX: doesn't plot line when no exons around

makeGeneTrack <- function(d.range, chr, collapseTx, exon.ref) {
    
    ## get the exons and gene IDs found in plotting range
    
    importRange <- GRanges(seqnames = chr, ranges = IRanges(start = d.range[1], end = d.range[2]))
    tx <- transcriptsByOverlaps(TxDb.Mmusculus.UCSC.mm10.knownGene, ranges = importRange, columns = c("TXNAME", "GENEID"))
    
    ## load in tx.sum df, y coord cap for y axis, and relative height of gene track    
    ## if collapseTx = TRUE, find the complete gene ranges of those found in plotting range
    if (collapseTx){
        gns <- unique(as.character(tx$GENEID))
        ex.df <- exon.ref[which(exon.ref$txname %in% gns),]
        tx.ref <- readRDS("./ref/mm10_TxDbKnownGene_TxCollapsedByGeneID_noScaffoldChrGeneIDonly_withYCoord_181118.rds")
        tx.sum <- tx.ref[which(tx.ref$geneid %in% gns),]
    } else {
        tx <- unique(tx$TXNAME)
        ex.df <- exon.ref[which(exon.ref$txname %in% tx),]
        ex.df <- ex.df[complete.cases(ex.df$geneid),]
        tx.sum <- suppressMessages(select(TxDb.Mmusculus.UCSC.mm10.knownGene, keys = tx, columns = c("GENEID", "TXSTART", "TXEND", "TXSTRAND"), keytype = "TXNAME"))
        colnames(tx.sum) <- tolower(colnames(tx.sum))
        tx.sum$genename <- ex.df$genename[match(tx.sum$txname, ex.df$txname)]
        tx.sum$ycoord <- ex.df$ycoord[match(tx.sum$txname, ex.df$txname)]
    }
 

    
    if (collapseTx){
        ycap <- 2
    } else {
        ycap <- 10
    }
    
    if(nrow(ex.df) == 0){
        gene.track <- plot_ly(x = d.range[1]:d.range[2], y = NA, type = "scatter")
        return(list(plot = list(gene.track), cap = ycap + 0.75))
    } else {
        
    

        #print(head(ex.df))
        #print(collapseTx)

            
        maxy <- max(ex.df$ycoord)
        #print(maxy)
        if (maxy > ycap) {
            scaleFactor <- ycap/maxy
        } else {
            scaleFactor <- 1
        }
        
        
        ## re-calculate ycoords according to scale factor
        
        ex.df$ycoord <- ex.df$ycoord*scaleFactor
        tx.sum$ycoord <- tx.sum$ycoord*scaleFactor
        
        
        ## calculate coords/direction for arrows and lines
        
        arrowBins <- 30
        points.coords <- seq(d.range[1], d.range[2], length.out = arrowBins)
        points.lst <- lapply(1:nrow(tx.sum), function(x){
            p.x <- points.coords[which(points.coords > tx.sum$txstart[x] & points.coords < tx.sum$txend[x])]
            data.frame(xcoord = p.x, 
                       txname = rep(tx.sum$txname[x], length(p.x)), 
                       ycoord = rep(tx.sum$ycoord[x], length(p.x)),
                       strand = rep(tx.sum$txstrand[x], length(p.x)),
                       gene = rep(tx.sum$genename[x], length(p.x)),
                       stringsAsFactors = FALSE)
            
        })
        pts.df <- do.call(rbind, points.lst)
        m <- tx.sum[,c("txstart", "txend")]
        m$filler <- NA
        g.lines <- as.vector(t(m))
        lines.df <- data.frame(xcoord = g.lines, 
                               ycoord = rep(tx.sum$ycoord, each = 3)
                               )
        
        lookup <- c("CDS", "utr|ncRNA")
        rect.lst <- lapply(lookup, function(x){
            df <- ex.df[grep(x, ex.df$regiontype),]
            if(nrow(df) > 0){
                m <- df[,c("start", "end")]
                m$filler <- NA
                g.lines <- as.vector(t(m))
                data.frame(xcoord = g.lines, ycoord = rep(df$ycoord, each = 3), gene = rep(df$genename, each = 3))
            } else {
                return(data.frame())
            }
    
        })
        names(rect.lst) <- lookup
    
        # ## make text annotation lists for plotly
        # 
        # anno <- list(
        #     xref = "x",
        #     yref = "y",
        #     showarrow = FALSE
        # )
        # 
        # annos <- list()
        # 
        # for (i in 1:nrow(tx.sum)) {
        #     
        #     str <- tx.sum$strand[i]
        #     
        #     if (str == "-"){
        #         anno[["x"]] <- tx.sum$end[i]
        #         anno[["xanchor"]] <- "left"
        #     } else {
        #         anno[["x"]] <- tx.sum$start[i]
        #         anno[["xanchor"]] <- "right"
        #     }
        #     anno[["text"]] <- paste0("<i>", tx.sum$gene.name[i], "</i>")
        #     anno[["y"]] <- tx.sum$ycoord[i]
        #     annos <- c(annos, list(anno))
        #     
        #     # if (str == "-"){
        #     #     anno.df$xcoord[i] <- tx.sum$end[i]
        #     # } else {
        #     #     anno.df$xcoord[i] <- tx.sum$start[i]
        #     # }
        #     
        # }
        
        
        ## calculate plotting coords for exons and make shape annotation lists for plotly
        
        # exHt <- scaleFactor * 0.70
        
        # rect <- list(
        #     type = "rect",
        #     opacity = 1,
        #     line = list(color = "black", width = 0.5),
        #     xref = "x",
        #     yref = "y"
        # )
        
        # rects <- list()
        # for (i in 1:nrow(ex.df)) {
        #     rect[["x0"]] <- ex.df$start[i]
        #     rect[["x1"]] <- ex.df$end[i]
        #     
        #     type <- ex.df$regiontype[i]
        #     
        #     if (type == "CDS") {
        #         ex.scale <- 1
        #         rect[["fillcolor"]] <- "black"
        #     } else {
        #         ex.scale <- 0.5
        #         rect[["fillcolor"]] <- "grey"
        #     }
        #     rect[c("y0")] <- ex.df$ycoord[i] - ex.scale * exHt/2
        #     rect[c("y1")] <- ex.df$ycoord[i] + ex.scale * exHt/2
        #     
        #     rects <- c(rects, list(rect))
        #     
        #     # gene.track <- gene.track %>% add_lines(
        #     #     x = c(ex.df$start[i], ex.df$end[i], ex.df$end[i], ex.df$start[i]),
        #     #     y = rep(c(ex.df$ycoord[i] - ex.scale * height/2, ex.df$ycoord[i] + ex.scale * height/2), 2),
        #     #     type = 'scatter', mode = "lines", 
        #     #     fill = 'toself',
        #     #     fillcolor = 'black',
        #     #     line = list(color = "black", width = 2)
        #     # )
        #     
        # }
        
        ## make plotly gene track plot (will add exon annotations later in layout)
        
        gene.track <- plot_ly() %>%
            add_trace(x = pts.df$xcoord, y = pts.df$ycoord, 
                      type = "scatter", mode = "markers",
                      text = pts.df$gene,
                      hoverinfo = "text", 
                      symbol = pts.df$strand, 
                      color = I("black"), 
                      fill = I("white"),
                      symbols = c("+" = "triangle-right-open", "-" = "triangle-left-open", "*" = "diamond-open"),
                      marker = list(size = 5)) %>%
            add_lines(x = lines.df$xcoord, y = lines.df$ycoord, 
                      color = I("black"), 
                      line = list(width = 0.5))
            
        
        if (nrow(rect.lst$CDS) > 0){
            gene.track <- gene.track %>%
                add_lines(x = rect.lst$CDS$xcoord, 
                          y = rect.lst$CDS$ycoord,
                          text = rect.lst$CDS$gene,
                          hoverinfo = "text",
                          color = I("black"), 
                          line = list(width = 10, simplify = FALSE))
        }
        
        if (nrow(rect.lst$utr) > 0){
            gene.track <- gene.track %>%
                add_lines(x = rect.lst$utr$xcoord, 
                          y = rect.lst$utr$ycoord, 
                          text = rect.lst$utr$gene,
                          hoverinfo = "text",
                          color = I("#7f7f7f"), 
                          line = list(width = 6, simplify = FALSE))
        }
        
        anno.df <- tx.sum
        anno.df$xcoord <- NA
        anno.lst <- lapply(c("+", "-", "*"), function(x) {
            df <- anno.df[which(anno.df$txstrand %in% x),]
            df$xcoord <- switch(x, 
                                `+` = df$txstart,
                                `-` = df$txend,
                                `*` = df$txstart)
            return(df)
        })
        names(anno.lst) <- c("middle left", "middle right", "middle left")
        for (i in 1:length(anno.lst)){
            df <- anno.lst[[i]]
            if (nrow(df) > 0){
                gene.track <- gene.track %>% add_trace(
                    x = df$xcoord,
                    y = df$ycoord,
                    text = paste0(" <i>", df$genename, "</i> "),
                    hoverinfo = "none",
                    type = 'scatter', mode = "text",
                    textposition = names(anno.lst)[i]
                )
            }
        }
        
        # return(list(plot = list(gene.track), rect = rects, anno = annos, ycap = ycap, height = plotHeight))
        return(list(plot = list(gene.track), cap = ycap + 0.75))
        
        # ## old code for making line annotations; instead made lines a trace
        # line <- list(
        #     type = "line",
        #     line = list(color = "black"),
        #     xref = "x",
        #     yref = "y"
        # )
        # 
        # lines <- list()
        # for (i in 1:nrow(tx.sum)) {
        #     line[["x0"]] <- tx.sum$start[i]
        #     line[["x1"]] <- tx.sum$end[i]
        #     line[c("y0", "y1")] <- tx.sum$ycoord[i]
        #     lines <- c(lines, list(line))
        #     
        #     # gene.track <- gene.track %>% add_lines(
        #     #     x = c(tx.sum$end[i], tx.sum$start[i]),
        #     #     y = rep(tx.sum$ycoord[i], 2),
        #     #     type = 'scatter', mode = "lines",
        #     #     line = list(color = "black", width = 0.5)
        #     # )
        # }
        
        
        # shapes <- c(lines, rects)
    }
}