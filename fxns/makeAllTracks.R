

# if plot just rendered, d == NULL
# if user double clicks for autozoom, then d$xaxis.autorange and d$yaxis.autorange == TRUE
# if user resizes page, then d$width is pixels of plot width
# if user drags with dragmode pan, then d$xaxis.range[0], d$xaxis.range[1] shows new start and end




makeAllTracks <- function(d.range, p.range, chr, collapseTx, exon.ref, st, type) {
        print(paste0("plotting using ", paste0(d.range, collapse = " "), " data range"))
        print(paste0("plotting using ", paste0(p.range, collapse = " "), " plotting range"))
        #d.range <- floor(d.range)
        #p.range <- flor(p.range)
    
        # make gene track ##############################################################
        gene.track <- NULL
        g <- makeGeneTrack(d.range = d.range, chr = chr, collapseTx, exon.ref)
        gene.track <- g$plot
        
        
        # get plotting elements for pileups and annotation ####################################
        
        fileLst <- split(st$bwFile, st$experiment)
        
        data.tracks <- c()
        yrange.lst <- c()
        anno.track <- plot_ly(type = "scatter", mode = "lines")
        
        for (i in 1:length(fileLst)) {
            expName <- names(fileLst)[i]
            col <- st$colors[match(expName, st$experiment)]
            lst <- fileLst[[i]]
            
            d.lst <- list()
            ymax.lst <- c()
            for (f in lst){
                p <- makePileupTrack(file = f, d.range = d.range,  p.range = p.range, col = col, chr = chr, type = type)
                d.lst <- c(d.lst, list(p$plot))
                ymax.lst <- c(ymax.lst, p$ymax)
            }
            
            ymax <- max(ymax.lst)
            toAdd <- rep(ymax, length(lst))
            names(toAdd) <- sapply(lst, function(x){
                info <- st[match(x, st$bwFile),]
                name <- paste0(info$sample, "\n", info$replicate)
                
                #gsub("_", "\n", name, fixed = TRUE)
            })
            yrange.lst <- c(yrange.lst, toAdd)
            data.tracks <- c(data.tracks, d.lst)
            
            dt.num <- nrow(st)
            annoFile <- unique(st$annoFile[match(expName, st$experiment)])
            
            anno.df <- makeAnnotationTrack(annoFile = annoFile, chr = chr, d.range = d.range,
                                           y.index = -i, col = col, yaxis = dt.num + 2)
            anno.track <- anno.track %>% add_lines(x = anno.df$xcoord, y = anno.df$ycoord, 
                                                   text = anno.df$text, hoverinfo = "text",
                                                   line = list(width = 15, simplify = FALSE, color = col))
        }
        
        print(paste0("this is y ranges: ", paste0(yrange.lst, collapse = ", ")))
        
        select.track <- plot_ly(x = seq(d.range[1], d.range[2]), y = 1, 
                                type = "scatter", mode = "markers",
                                marker = list(color = "black", size = 1, symbol = "square")) %>% toWebGL()
        
        ## generate subplots ################################
        
        
        
        s <- c(gene.track, data.tracks, list(anno.track), list(select.track))
        s.heights <- c(0.2, rep((0.6)/length(data.tracks), length(data.tracks)), 0.1, 0.1)
        
        sl <- subplot(s,  nrows = length(s), shareX = TRUE, heights = s.heights) %>%
            layout(dragmode = "pan", showlegend = FALSE,
                   xaxis = list(title = "genomic coordinates (bp)", range = c(p.range[1], p.range[2])),
                   yaxis = list(title = "y1", type = "linear", fixedrange = TRUE, range = c(0, g$cap), visible = FALSE))
        
        for (i in 1:length(data.tracks)){
            cap <- yrange.lst[i]
            name <- names(yrange.lst)[i]
            l <- list(list(title = name, type = "linear", fixedrange = TRUE, range = c(0, cap)))
            names(l) <- paste0("yaxis", i+1)
            arg <- c(list(sl), l)
            sl <- do.call("layout", arg)
        }

        al <- list(list(title = "peak\nregions", 
                        type = "linear", 
                        fixedrange = TRUE, 
                        showline = FALSE,
                        showticklabels = FALSE,
                        ticks = "")
                   )
        names(al) <- paste0("yaxis", length(s)-1)
        arg <- c(list(sl), al)
        sl <- do.call("layout", arg)
        
        sel.l <- list(list(title = "select\nhere", 
                           fixedrange = TRUE, 
                           showline = FALSE,
                           showticklabels = FALSE,
                           ticks = "")
                      )
        names(sel.l) <- paste0("yaxis", length(s))
        arg <- c(list(sl), sel.l)
        sl <- do.call("layout", arg)
        
        sl %>% config(collaborate = FALSE, displaylogo = FALSE, editable = TRUE,
                     modeBarButtonsToRemove= c("sendDataToCloud", 
                                               "lasso2d", 
                                               "autoScale2d", 
                                               "resetScale2d", 
                                               "hoverClosestCartesian", 
                                               "hoverCompareCartesian",
                                               "toggleSpikelines"))
}