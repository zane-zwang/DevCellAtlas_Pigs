# Purpose: plot relative expression helper.

## PlotRE
PlotRE_z <- function(ExpressionSet,
                   Groups     = NULL,
                   modules    = NULL,
                   legendName = "age",
                   xlab = "Ontogeny",
                   ylab = "Relative Expression Level",
                   main = "",
                   y.ticks = 10,
                   adjust.range = TRUE,
                   alpha = 0.008, ...)
{
        
        ExpressionSet <- as.data.frame(ExpressionSet)
        is.ExpressionSet(ExpressionSet)
        
        stage <- expr <- age <- NULL
        
        if(is.null(Groups))
                stop("Your Groups list does not store any items.", call. = FALSE)
        
        ### getting the PS names available in the given expression set
        age_names <- as.character(names(table(ExpressionSet[ , 1])))
        
        # Require every group element to have an age value.
        if(!all(unlist(Groups) %in% as.numeric(age_names)))
                stop("There are items in your Group elements that are not available in the age column of your ExpressionSet.", call. = FALSE)
        
        if (length(Groups) > 2)
                stop("Please specify at maximum 2 groups that shall be compared.", call. = FALSE)
        
        ### getting the PS names available in the given expression set
        nPS <- length(age_names)
        nCols <- dim(ExpressionSet)[2]
        ### define and label the REmatrix that holds the rel. exp. profiles
        ### for the available PS
        MeanValsMatrix <- matrix(NA_real_,nPS,nCols-2)
        rownames(MeanValsMatrix) <- age_names
        colnames(MeanValsMatrix) <- names(ExpressionSet)[3:nCols]
        
        MeanValsMatrix <- age.apply(ExpressionSet, RE)
        mean.age <- data.frame(age = age_names, MeanValsMatrix, stringsAsFactors = FALSE)
        mMatrix <- tibble::as_tibble(reshape2::melt(mean.age, id.vars = "age"))
        colnames(mMatrix)[2:3] <- c("stage", "expr")
        
        if (length(Groups) == 1) {
                
                p <- ggplot2::ggplot(mMatrix, ggplot2::aes( factor(stage, levels = unique(stage)), expr, group = age, fill = factor(age, levels = age_names))) + 
                        ggplot2::geom_line(ggplot2::aes(color = factor(age, levels = age_names)), size = 2) +
                        ggplot2::labs(x = xlab, y = ylab, title = main, colour = legendName) +
                        ggplot2::theme_classic() +
                        ggplot2::theme(
                                aspect.ratio = 0.55,
                                title            = ggplot2::element_text(size = 12),
                                legend.title     = ggplot2::element_text(size = 12),
                                legend.text      = ggplot2::element_text(size = 12),
                                axis.title       = ggplot2::element_text(size = 12),
                                axis.text.y      = ggplot2::element_text(size = 12),
                                axis.text.x      = ggplot2::element_text(size = 12),
                                panel.background = ggplot2::element_blank(),
                                strip.text.x     = ggplot2::element_text(
                                        size           = 12,
                                        colour         = "black"
                                )
                        ) +
                        ggplot2::scale_y_continuous(breaks = scales::pretty_breaks(n = y.ticks)) + 
                        ggplot2::scale_colour_manual(values = custom.myTAI.cols(nrow(mMatrix))) 
                
                
                if (!is.null(modules)) {
                        p <- p + ggplot2::geom_rect(data = mMatrix,ggplot2::aes(
                                xmin = modules[[2]][1],
                                xmax = modules[[2]][length(modules[[2]])],
                                ymin = min(MeanValsMatrix) - (min(MeanValsMatrix) / 50),
                                ymax = Inf), fill = "#4d004b", alpha = alpha)  
                }
                p <- p + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 1,hjust = 1))
                return(p)
        }
        
        if (length(Groups) == 2) {
                
                mMatrixGroup1 <- dplyr::filter(mMatrix, age %in% Groups[[1]])
                mMatrixGroup2 <- dplyr::filter(mMatrix, age %in% Groups[[2]])
                
                p1 <- ggplot2::ggplot(mMatrixGroup1, ggplot2::aes( factor(stage, levels = unique(stage)), expr, group = age, fill = factor(age, levels = age_names[Groups[[1]]]))) + 
                        ggplot2::geom_line(ggplot2::aes(color = factor(age, levels = age_names[Groups[[1]]])), size = 2) +
                        ggplot2::labs(x = xlab, y = ylab, title = main, colour = legendName) +
                        ggplot2::theme_classic() +
                        ggplot2::theme(
                                aspect.ratio = 0.55,
                                title            = ggplot2::element_text(size = 12),
                                legend.title     = ggplot2::element_text(size = 12),
                                legend.text      = ggplot2::element_text(size = 12),
                                axis.title       = ggplot2::element_text(size = 12),
                                axis.text.y      = ggplot2::element_text(size = 12),
                                axis.text.x      = ggplot2::element_text(size = 12),
                                panel.background = ggplot2::element_blank(),
                                strip.text.x     = ggplot2::element_text(
                                        size           = 12,
                                        colour         = "black"
                                )
                        ) +
                        ggplot2::scale_colour_manual(values = custom.myTAI.cols(nrow(mMatrix))[Groups[[1]]]) +
                        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = -90, hjust = 0))
                if (!adjust.range) {
                        
                        p1 <- p1 + ggplot2::scale_y_continuous(breaks = scales::pretty_breaks(n = y.ticks))
                }
                
                if (!is.null(modules)) {
                        p1 <- p1 + ggplot2::geom_rect(data = mMatrixGroup1, ggplot2::aes(
                                xmin = modules[[2]][1],
                                xmax = modules[[2]][length(modules[[2]])],
                                ymin = min(MeanValsMatrix),
                                ymax = Inf), fill = "#4d004b", alpha = alpha)  
                }
                
                p2 <- ggplot2::ggplot(mMatrixGroup2, ggplot2::aes( factor(stage, levels = unique(stage)), expr, group = age, fill = factor(age, levels = age_names[Groups[[2]]]))) + 
                        ggplot2::geom_line(ggplot2::aes(color = factor(age, levels = age_names[Groups[[2]]])), size = 2) +
                        ggplot2::labs(x = xlab, y = ylab, title = main, colour = legendName) +
                        ggplot2::theme_minimal() +
                        ggplot2::theme(
                                aspect.ratio = 0.55,
                                title            = ggplot2::element_text(size = 12),
                                legend.title     = ggplot2::element_text(size = 12),
                                legend.text      = ggplot2::element_text(size = 12),
                                axis.title       = ggplot2::element_text(size = 12),
                                axis.text.y      = ggplot2::element_text(size = 12),
                                axis.text.x      = ggplot2::element_text(size = 12),
                                panel.background = ggplot2::element_blank(),
                                strip.text.x     = ggplot2::element_text(
                                        size           = 12,
                                        colour         = "black"
                                )
                        ) + 
                        ggplot2::scale_colour_manual(values = custom.myTAI.cols(nrow(mMatrix))[Groups[[2]]]) +
                        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = -90, hjust = 0))
                
                if (!adjust.range) {
                        
                        p2 <- p2 + ggplot2::scale_y_continuous(breaks = scales::pretty_breaks(n = y.ticks)) 
                }
                
                
                if (adjust.range){
                        p1 <- p1 + ggplot2::scale_y_continuous(limits = c(min(MeanValsMatrix), max(MeanValsMatrix)), breaks = scales::pretty_breaks(n = y.ticks))
                        
                        p2 <- p2 + ggplot2::scale_y_continuous(limits = c(min(MeanValsMatrix), max(MeanValsMatrix)), breaks = scales::pretty_breaks(n = y.ticks))    
                }
                
                if (!is.null(modules)) {
                        p2 <- p2 + ggplot2::geom_rect(data = mMatrixGroup2,ggplot2::aes(
                                xmin = modules[[2]][1],
                                xmax = modules[[2]][length(modules[[2]])],
                                ymin = min(MeanValsMatrix),
                                ymax = Inf), fill = "#4d004b", alpha = alpha)  
                }
                p1 <- p1 + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 1,hjust = 1))
                p2 <- p2 + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 1,hjust = 1))
                return(ggpubr::ggarrange(p1, p2, ncol = 2))
        }
}
