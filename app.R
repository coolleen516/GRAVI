# Load libraries ##############################################################

suppressMessages(library(shiny))
suppressMessages(library(plotly))
suppressMessages(library(dplyr))
suppressMessages(library(rtracklayer))
suppressMessages(library(RColorBrewer))
suppressMessages(library(TxDb.Mmusculus.UCSC.mm10.knownGene))
suppressMessages(library(zoo)) # to calculate AUC
options(stringsAsFactors = FALSE)


# Load in global objects ######################################################

# annotation reference for matching gene IDs and gene names
anno.ref <- readRDS("./ref/mm10_TxDbKnownGene_geneIDtoNameDescr_170526.rds")
gns <- as.character(anno.ref$gene.name)

# chromosome data for generating ideograms
chrom.data <- read.delim("./ref/UCSC_mm10_cytoBandIdeo_181124.txt")
colnames(chrom.data)[1] <- "chrom"


# Load in functions ###########################################################

fxns <- list.files("./fxns", full.names = TRUE)

for (fxn in fxns){
    source(fxn)
}

# Auto-detect datasets ########################################################

sampletable <- makeSampleTable()
choices <- unique(sampletable$experiment)
names(choices) <- gsub("_", " ", choices)

# Load in colors options for plotting #########################################

colors <- c(brewer.pal(9, "Set1")[1:5], rev(brewer.pal(8, "Dark2")))
colornames <- c("Fire Engine Red", "Curious Blue", "Fruit Salad", 
                "Violet Blue", "Dark Orange", "Dim Gray", "Hot Toddy", 
                "Gamboge", "Christi", "Deep Cerise", "Chetwode Blue", 
                "Tawny", "Elf Green")
color.key <- data.frame(colors = colors, colornames = colornames, exp = NA)
appCSS <- paste0(".colSelect .selectize-control .option:nth-child(", 
                 1:length(colors), ") { background-color: ", colors,
                 "; color: white }", collapse = " ")


# Define UI for application ###################################################

ui <- fluidPage(
    
    # load in CSS style for color option inputs
    tags$head(
        tags$style(appCSS)
    ),
    
    # application title
    headerPanel(
        list( img(src="Gravy_Boat-icon.png", height = "50px", width = "50px"), 
              strong("GRAVI", style = "font-size:30px"),
              strong("(v.1.0)", style = "font-size:18px")
        ),
        windowTitle="GRAVI"
    ),
    div(HTML("<strong>G</strong>enomic <strong>R</strong>egulation 
             <strong>A</strong>nalyzed <strong>V</strong>isually and 
             <strong>I</strong>nteractively")),
    br(),
   
    # this script unfocuses the buttons when pressed
    tags$script(HTML("$(document).ready(function() {
                     $('.btn').on('click', function(){$(this).blur()});
                     })")
                   ),
    

    fluidRow(
        
        # side panel for dataset and plotting options
        column(width = 3,
               wellPanel(
                   radioButtons("collapseTx", "Gene display: ", 
                                c("Expanded" = FALSE, "Collapsed" = TRUE), 
                                selected = TRUE)
               ),
               wellPanel(
                   
                   selectInput("dataset", label = "Pick Seq Data",
                               choices = as.list(choices),
                               selected = choices[1]
                               ),
                   selectInput("sample", label = "Pick condition",
                               choices = as.list(choices),
                               selected = choices[1]
                               ),
                   selectInput("replicate", label = "Pick replicate",
                               choices = as.list(choices),
                               selected = choices[1]
                               ),
                   tags$div(
                       tags$style(
                           HTML('#add{background-color:#4c83b6; color:white}
                                 #clear{background-color:lightgrey'
                           )
                       )
                   ),
                   actionButton("add", "Add samples"),
                   actionButton("clear", "Clear samples")
               ),
               wellPanel(
                   selectizeInput("genename", label = "Gene name", 
                                  choices = NULL),
                   tags$div(
                       tags$style(
                           HTML('#submit{background-color:#4c83b6; color:white}'
                           )
                       )
                   ),
                   actionButton("submit", "Go to gene"),
                   div(id = "errMsgPlace1")                   
               ),
               wellPanel(
                   tags$div(
                       tags$style(
                           HTML('#replot{background-color:#4c83b6; color:white}'
                           )
                       )
                   ),
                   actionButton("replot", label = "  Re-plot", 
                                icon = icon("refresh"), width = "100%"),
                   div(id = "errMsgPlace3")
               ),
               uiOutput("traceColors")
      ),
      
      
      # Tabset for sampletable display and plots
      column( width = 9,
          tabsetPanel( id = "inTabset",
              tabPanel("Sample table", value = "tabSt",
                       fluidRow(
                           column(width = 8,
                               tableOutput("table")
                           )
                           
                       )
              ),
              tabPanel("Plot", value = "tabPlot",
                       column(width = 8,
                              br(),
                              plotlyOutput("tracks", height = "800px"),
                              plotlyOutput("ideogram", height = "200px")
                              ),
                       column(width = 4,
                              br(),
                              br(),
                              div(id = "errMsgPlace2"),
                              plotlyOutput("barplots", height = "600px")
                              )
              ),         
              tabPanel("Extra", value = "tabExtra",
                       verbatimTextOutput("zoom"), 
                       verbatimTextOutput("selected"))
              )

      )

   )

)

# Define server logic #########################################################

server <- function(input, output, session) {
    
    
    ## automatically load datasets ############################################
    
    observe({
        sub.st <- sampletable[which(sampletable$experiment %in% input$dataset),]
        samples <- unique(sub.st$sample)
        names(samples) <- samples
        updateSelectInput(session, "sample",
                          choices = as.list(samples),
                          selected = samples[1]
        )
    })
    
    observe({
        sub.st <- sampletable[which(sampletable$experiment %in% input$dataset & 
                                    sampletable$sample %in% input$sample),]
        reps <- c(unique(sub.st$replicate), "ALL")
        names(reps) <- reps
        updateSelectInput(session, "replicate",
                          choices = as.list(reps),
                          selected = 1
        )
        
    })
    
    observe({
        updateSelectizeInput(session, 'genename', choices = gns, server = TRUE)
    })
    
    # load in sampletable #####################################################
    ## includes experiment, sample, replicate, paths for bigwig (bwFile) and 
    ## peak files (annoFile), and assigned colors
    
    
    # reactive values: st = sampletable, colKey = color key, uiKey = UI key for
    # color options
    samples <- reactiveValues(st = NULL, colKey = NULL, uiKey = NULL)
    
    
    # when add is pushed, samples and info are added row by row, and colors are
    # are by color order
    observeEvent(input$add, {
        
        cat(file=stderr(), "___Step1: Loading samples___", "\n")

        sample.df <- sampletable[which(sampletable$experiment %in% 
                                       input$dataset & sampletable$sample %in% 
                                       input$sample),]

        if (input$replicate != "ALL") {
            sample.df <- sample.df[which(sample.df$replicate %in% 
                                        input$replicate),]
        }
        
        # assign colors in color key
        st <- samples$st
        if (is.null(st)){
            newExp <- unique(sample.df$experiment)
            color.key[1,"exp"] <- newExp
            samples$colKey <- color.key # sets initial color key if first

        } else {
            currExp <- unique(st$experiment)
            newExp <- input$dataset
            cat(file=stderr(), paste0("currExp: ", currExp), "\n")
            if (!(newExp %in% currExp)){
                colKey <- samples$colKey
                
                # assign color to whatever the first NA exp is
                colKey[which(is.na(colKey$exp))[1], "exp"] <- newExp
                samples$colKey <- colKey
            }
        }
        
        curr.key <- samples$colKey[complete.cases(samples$colKey),]
        curr.key$ids <- paste0("colSelect", 1:nrow(curr.key))
        samples$uiKey <- curr.key
        
        cat(file=stderr(), paste0("newExp: ", newExp), "\n")
        cat(file=stderr(), paste0("colKey: ", "\n"))
        cat(file=stderr(), paste0(samples$colKey$colornames, "\n"))
        
        # append samples to samples$st reactive value
        sample.df$colors <- samples$colKey$colors[match(sample.df$experiment, 
                                                        samples$colKey$exp)]
        samples$st <- isolate ({rbind (samples$st, sample.df)})
        
    })

    
    # generate new UI for each dataset ########################################
    output$traceColors <- renderUI({
        
        req(samples$uiKey)
        isolate({ 
            
            lst <- list()
            cat(file=stderr(), "___Step2: Making UI___", "\n")
            cols <- setNames(color.key$colors, color.key$colornames)
            
            for (i in 1:nrow(samples$uiKey)){
                
                exp <- samples$uiKey$exp[i]
                col <- samples$uiKey$colors[i]
                id <- samples$uiKey$ids[i]
                
                lst[[i]] <- div(selectInput(inputId = id, 
                                            label = paste0("Colors for ", 
                                                           exp, " traces"),
                                            choices = cols, selected = col), 
                                class = "colSelect"
                                )
            }
            
            tagList(lst)
            
        })
        
    })
    
    # if new inputs, generate st with updated colors ##########################
    observe({
        
        req(samples$st)
        ui.key <- samples$uiKey
        st <- samples$st
        
        cat(file=stderr(), "___Step3: Updating colors with input for ids___", "\n")
        
        for (i in ui.key$ids){
            
            req(input[[i]])
            cat(file=stderr(), paste0( i, "\n"))
            ui.key$colors[which(ui.key$ids %in% i)] <- input[[i]]
            updateSelectInput(session, inputId = i, selected = input[[i]])
            
        }
        
        st$colors <- ui.key$colors[match(st$experiment, ui.key$exp)]
        samples$st <- st
        cat(file=stderr(), "Updated st: ", "\n")
        cat(file=stderr(),  st$colors, "\n")

    })

    
    # clear samples by making samples$df NULL #################################
    observeEvent(input$clear, {
        
        samples$st <- NULL
        cat(file=stderr(), "___Clear samples___", "\n")
        
    })
    
    
    # render sampletable ######################################################
    output$table <- renderTable({
        
        input$add
        
        if(!is.null(samples$st) > 0){
            df <- samples$st
            df[,c("experiment", "sample", "replicate")]
        } else {
            data.frame(experiment = "", sample = "", replicate = "")
        }
    })
    
    # find range by gene ######################################################
    
    # reactive values for index/location of bp ranges and gene-related items
    ind <- reactiveValues(p.range = NULL, d.range = NULL, reset = NULL)
    gene <- reactiveValues(generange = NULL, chr = NULL, exon.ref = NULL)
    
    
    # load gene-related parameters
    observeEvent(input$submit, {
        
        cat(file=stderr(), "___Step4: Plotting based on gene location___", "\n")
        
        # error message if no samples are loaded
        if (is.null(samples$st)){
            insertUI(
                selector = "#errMsgPlace1",
                where = "afterEnd",
                ui = p("Please add samples first.", style = "font-size:12px; color:red", id = "errMsg1")
            )
            req(FALSE)
        } else {
            removeUI(
                selector = "#errMsg1"
            )
        }
        
        # switch to 'plot' tab
        updateTabsetPanel(session, "inTabset",
                          selected = "tabPlot")
        
        
        # load in exon ref
        collapseTx <- input$collapseTx
        exon.ref <- switch(collapseTx,
                           `FALSE` = readRDS( file="./ref/mm10_TxDbKnownGene_ExonTypes_noScaffoldChrGeneIDonly_withYCoord_181117.rds"),
                           `TRUE` = readRDS( file ="./ref/mm10_TxDbKnownGene_ExonCollapsedByGeneID_noScaffoldChrGeneIDonly_withYCoord_181117.rds")
        )
        
        gene$exon.ref <- exon.ref 
        
        # get gene coordinates
        geneID <- anno.ref$gene.id[which(anno.ref$gene.name %in% input$genename)]
        gn.ex.df <- exon.ref[which(exon.ref$geneid %in% geneID),]
        chr <- unique(as.character(gn.ex.df$txchrom))
        gene$chr <- chr #*************
        gn.start <- min(gn.ex.df$start)
        gn.end <- max(gn.ex.df$end)
        
        # p.range is plotting range
        gene.w <- gn.end - gn.start
        gene.range <- c(gn.start - floor(0.1*gene.w), gn.end + floor(0.1*gene.w))
        gene$generange <- gene.range
        p.range <- gene.range
        cat(file=stderr(), "loading initial plotting range: ", paste0(p.range, collapse = ", "), "\n")
        ind$p.range <- p.range
        
        # d.range is data range
        w <- diff(p.range)
        d.range <- c(p.range[1] - w, p.range[2] + w)
        cat(file=stderr(), paste0("loading initial data range: ", paste0(d.range, collapse = ", ")), "\n")
        ind$d.range <- d.range
        
        # trigger to plot by gene
        ind$reset <- TRUE
        
    })
    
    ## generate tracks ########################################################
    ## by modes: triggered by submit, out of range data, or manual replotting
    
    # this function is here to plot 
    # (1) when submit is pressed (input$reset = TRUE)
    # (2) auto plot if when dragging/zooming out of data range
    # Note: (2) requires that it is always aware of d.range and p.range
    
    observe({
        
        cat(file=stderr(), "___Step5: Plot gene, new boundaries, or load new ranges___", "\n")
        
        req(ind$d.range)
        req(ind$p.range)
        
        st <- isolate({samples$st})
        collapseTx <- isolate({input$collapseTx})
        exon.ref <- isolate({gene$exon.ref})
        chr <- isolate({gene$chr})
        test <- isolate({ind$reset})
        
        # load in event data of plot (i.e. plot coords)
        d <- event_data("plotly_relayout", source = "A")
        
        if (test == TRUE){
            
            cat(file=stderr(), "plotting gene location", "\n")
            
            p.range <- ind$p.range
            d.range <- ind$d.range
            ind$reset <- FALSE
            
            output$tracks <- renderPlotly({
                
                progress <- Progress$new(session, min=1, max=15)
                on.exit(progress$close())

                progress$set(message = 'Plotting',
                             detail = 'Doo doo doo...')
                
                makeAllTracks(d.range, p.range, chr, collapseTx, exon.ref, st,
                              type = "line")
                
            })
            
        } else {
            
            d.range <- isolate({ind$d.range})
            ind$reset <- FALSE
            
            # if d has plot coordinates, then only replot if out of d.range
            if (all(c("xaxis.range[0]", "xaxis.range[1]") %in% names(d))){
                
                p.range <- floor(c(d$`xaxis.range[0]`,d$`xaxis.range[1]`))
                w <- diff(p.range)
                cat(file=stderr(), paste0("detected new p.range: ", 
                             paste0(p.range, collapse = " ")), "\n")
                
                ## if new data outside of boundaries, then load and plot.
                ## otherwise, load in new coords and end observation
                if (ind$d.range[2] < p.range[2] | ind$d.range[1] > p.range[1]) {
                    
                    d.range <- c(p.range[1] - w, p.range[2] + w)
                    ind$d.range <- d.range
                    ind$p.range <- p.range

                    output$tracks <- renderPlotly({
                        
                        cat(file=stderr(), "plotting and loading new boundaries", "\n")
                        
                        progress <- Progress$new(session, min=1, max=15)
                        on.exit(progress$close())
                        
                        progress$set(message = 'Extending tracks',
                                     detail = 'Yay...')
                        
                        makeAllTracks(d.range, p.range, chr, collapseTx, 
                                      exon.ref, st, type = "line")
                        
                    })
                    
                } else {

                    ind$p.range <- p.range
                    cat(file=stderr(), paste0("loaded in new plotting coords, ", 
                                 paste0(p.range, collapse = " ")), "\n")
                    req(FALSE)
                    
                }
            
            # if d does not have x coords (i.e. if you just select something
            # else on modebar), then end observation      
            } else {
                
                cat(file=stderr(), "no change to p.range", "\n")
                req(FALSE)
                
            }
            
        }

    })
    
    
    # trigger manual replotting
    # used for things like:
    # (1) higher resolution of peaks when zoomed in
    # (2) changing color of peaks
    # (3) re-plotting when new samples are added
    # (4) re-plotting when gene display is changed
    
    observeEvent(input$replot, {
        
        # error message if samples are not loaded
        if (is.null(samples$st)){
            insertUI(
                selector = "#errMsgPlace3",
                where = "afterEnd",
                ui = p("Please add samples first.", 
                       style = "font-size:12px; color:red", 
                       id = "errMsg3")
            )
            req(FALSE)
        } else {
            removeUI(
                selector = "#errMsg3"
            )
        }
        
        st <- samples$st
        d <- event_data("plotly_relayout", source = "A")

        d.range <- ind$d.range
        p.range <- ind$p.range
        collapseTx <- input$collapseTx
        chr <- gene$chr
        exon.ref <- switch(collapseTx,
                           `FALSE` = readRDS( file="./ref/mm10_TxDbKnownGene_ExonTypes_noScaffoldChrGeneIDonly_withYCoord_181117.rds"),
                           `TRUE` = readRDS( file ="./ref/mm10_TxDbKnownGene_ExonCollapsedByGeneID_noScaffoldChrGeneIDonly_withYCoord_181117.rds")
        
        )
        gene$exon.ref <- exon.ref

        # if percentage of plotting range is less than 20% of data range, then reset data range
        w <- diff(p.range)
        if (w/diff(d.range) < 0.2){
            d.range <- c(p.range[1] - w, p.range[2] + w)
            ind$d.range <- d.range
        }
        
        output$tracks <- renderPlotly({
            
            cat(file=stderr(), "___Step6: Replotting___", "\n")
            cat(file=stderr(), paste0("p.range: , ", paste0(p.range, collapse = " ")), "\n")
            cat(file=stderr(), paste0("d.range: , ", paste0(d.range, collapse = " ")), "\n")
            
            progress <- Progress$new(session, min=1, max=15)
            on.exit(progress$close())
            
            progress$set(message = 'Re-plotting',
                         detail = 'Weeeeee...')
            
            makeAllTracks(d.range, p.range, chr, collapseTx, exon.ref, st, 
                          type = "line")
            
        })
        
    })
    
    
    # plot chromosome ideogram ################################################
    
    # first plot ideogram
    output$ideogram <- renderPlotly({
        
            cat(file=stderr(), "___Step7a: Ideogram___", "\n")
        
            req(gene$chr)
            chr <- gene$chr
            
            toPlot <- chrom.data[which(chrom.data$chrom %in% chr),]
            levs <- c("gneg", "gpos33", "gpos66", "gpos75", "gpos100")
            cols <- c(brewer.pal(length(levs), "Greys"))
            names(cols) <- levs
            max.x <- max(toPlot$chromEnd)
            
            p <- plot_ly(source = "C") %>% 
                    layout(title = chr,
                           margin = list(t = 50, r = 50),
                           yaxis = list(visible = FALSE, fixedrange = TRUE),
                           xaxis = list(visible = FALSE, fixedrange = TRUE, 
                                        range = c(0, max.x)),
                           showlegend = FALSE)
            
            for (i in 1:nrow(toPlot)) {
                
                x1 <- toPlot$chromStart[i]
                x2 <- toPlot$chromEnd[i]
                col <- cols[[toPlot$gieStain[i]]]
                p <- p %>% 
                        add_trace(x = c(x1,x1,x2,x2), y = c(0,1,1,0),
                                type = 'scatter',
                                mode = 'lines',
                                fill = 'tozeroy',
                                fillcolor = col,
                                line = list( color = 'black'),
                                text = toPlot$name[i],
                                hoveron = "fills",
                                hoverinfo = 'text'
                                )
                
            }
            
            p <- p %>% 
                add_trace(x = c(0, max.x), y = c(0,0), 
                         type = "scatter", 
                         mode = "lines", 
                         line = list(color = 'black'))
            
            p %>% 
                add_trace(x = c(ind$p.range[1],ind$p.range[1],
                                  ind$p.range[2],ind$p.range[2]),
                        y = c(-0.5,1.5,1.5,-0.5),
                        type = 'scatter',
                        mode = 'lines',
                        fill = 'toself',
                        fillcolor = "red",
                        line = list(color = 'red'),
                        text = "current",
                        hoveron = "fills",
                        hoverinfo = 'text') %>% 
                    config(displayModeBar = FALSE)
            
    })
    
    # then plot rectangle that highlights region
    observe({
        
        req(ind$p.range)
        p.range <- ind$p.range
        cat(file=stderr(), "___Step7b: Updating highlight of ideogram___", "\n")
        cat(file=stderr(), paste0("p.range start: ", p.range[1]), "\n")
        
        mid.pt <- ind$p.range[1] + round(diff(ind$p.range)/2)
        
        plotlyProxy("ideogram", session) %>%
            plotlyProxyInvoke("deleteTrace", list(-1))
        
        traceToAdd <- list(x = c(ind$p.range[1],ind$p.range[1],
                                 ind$p.range[2],ind$p.range[2]),
                           y = c(-0.5,1.5,1.5,-0.5),
                           type = 'scatter',
                           mode = 'lines',
                           fill = 'toself',
                           fillcolor = "red",
                           line = list(color = 'red'),
                           text = "current",
                           hoveron = "fills",
                           hoverinfo = 'text')
        
        layoutToAdd <- list(title = gene$chr,
                            margin = list(t = 50),
                            yaxis = list(visible = FALSE),
                            xaxis = list(visible = FALSE),
                            showlegend = FALSE)
        
        plotlyProxy("ideogram", session) %>%
            plotlyProxyInvoke("addTraces", traceToAdd) #%>%
 
    })
    
    # plot reactive barplots ##################################################
    
    output$barplots <- renderPlotly({
        
        d <- event_data("plotly_selected", source = "A")
        req(d)
        
        cat(file=stderr(), "___Step8: Barplots___", "\n")
        
        # error message if user does not select for designated region
        if(is.null(d) == T){
           return(NULL)
        } else if (class(d) == "list") {
            
               insertUI(
                   selector = "#errMsgPlace2",
                   where = "afterEnd",
                   ui = p("Please select an area from selection region.", style = "font-size:24px; color:grey", id = "errMsg2")
               )
               req(FALSE)
               
        } else {
            
           removeUI(
               selector = "#errMsg2"
           )
           
           xrange <- c(min(d$x), max(d$x))

           plt <- rbind()
           st <- isolate({samples$st})
           
           for (i in 1:nrow(st)){
               file <- st[i, "bwFile"]
               col <- st$colors[i]
               sampleID <- paste0(st$sample[i], "_", st$replicate[i])
               chr <- isolate({gene$chr})
               xy <- getBwScores(file = file, d.range = xrange, 
                                 maxbp = Inf, chr = chr)
               AUC <- sum(diff(xy$x)*rollmean(xy$y,2))
               df <- data.frame(x = sampleID, 
                                y = AUC, 
                                col = col, 
                                sample = st$sample[i], 
                                replicate = st$replicate[i],
                                experiment = st$experiment[i]
               )
               
               plt <- rbind(plt, df)
               
           }
           
           plt.sum <- plt %>% 
                        group_by(experiment, sample) %>% 
                            dplyr::summarise(mean = mean(y),
                                             sd = sd(y), 
                                             col = unique(col)
                                             )
           blnk <- plt
           blnk$y <- blnk$y * 1.1
           
           pltCols <- plt.sum$col
           names(pltCols) <- plt.sum$experiment
           pltCols <- unique(pltCols)
           plt.sum <<- plt.sum
           plt <<- plt
           blnk <<- blnk
           gg <- ggplot(data = plt.sum, aes(x = sample, 
                                            y = mean, 
                                            fill = experiment)) +
                 geom_bar(aes(label = mean), 
                          stat = "identity", 
                          color = "black") +
                 scale_fill_manual(values = pltCols) +
                 geom_errorbar(data = plt.sum, 
                               aes(ymin = mean - sd, ymax = mean + sd), 
                               width = 0.2, color = "black") +
                 geom_point(data = plt, aes( y = y, label = y), 
                            shape = 21, 
                            color = "black") +
                 geom_blank(data = blnk, aes(y = y)) +
                 scale_y_continuous(expand = c(0,0), limits = c(0, NA)) +
                 ylab("AUC") +
                 theme_classic()
           
           g <- ggplotly(gg, tooltip = c("label"), source = "B") %>%
                    layout(title = "Area under curve",
                          margin = list(t = 80),
                          yaxis = list(domain = c(0, 0.9)),
                          legend = list(orientation = "h", x = 0, y = -0.2))
           
           config(g, collaborate = FALSE, displaylogo = FALSE,
                  modeBarButtonsToRemove= c("pan2d",
                                            "select2d",
                                            "zoomIn2d",
                                            "zoomOut2d",
                                            "sendDataToCloud", 
                                            "lasso2d", 
                                            "autoScale2d", 
                                            "hoverClosestCartesian", 
                                            "hoverCompareCartesian",
                                            "toggleSpikelines"))
           
        }

    })
   
    
    # extra info about selections #############################################
    output$zoom <- renderPrint({
       d <- event_data("plotly_relayout", source = "A")
       if (is.null(d)) "Relayout (i.e., zoom) events appear here" else d
    })
    output$selected <- renderPrint({
       d <- event_data("plotly_selected", source = "A")
       if (is.null(d)) "Relayout (i.e., zoom) events appear here" else d
    })
   

   # stop session when browser closed #########################################
   session$onSessionEnded(function() {
       stopApp()
   })

   
}

# Run the application 
shinyApp(ui = ui, server = server, options = list(port=2345))
