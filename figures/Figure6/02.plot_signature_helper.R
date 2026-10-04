# Purpose: plot signature helper.

## PlotSignature_z_edit
PlotSignature_z <-
        function(ExpressionSet,
                 measure = "TAI",
                 TestStatistic = "FlatLineTest",
                 modules = NULL,
                 permutations = 1000,
                 lillie.test  = FALSE,
                 p.value = TRUE,
                 shaded.area  = FALSE,
                 custom.perm.matrix = NULL,
                 xlab = "Ontogeny",
                 ylab = "TAI",
                 main = "",
                 lwd = 2,
                 alpha = 0.1,
                 y.ticks = 10) {
                
        # [WARN] (Stefan) - Since all of TAI, TDI and TPI compute the same thing,
        # is there even a point of specifying the measure?
        # the ExpressionSet already specifies what type of index is being used
        
        if (!is.element(measure, c("TAI", "TDI", "TPI")))
                stop(
                        "Measure '",
                        measure,
                        "' is not available for this function. Please specify a measure supporting by this function.",
                        call. = FALSE
                )
        
        if (!is.element(TestStatistic, c("FlatLineTest", "ReductiveHourglassTest", "EarlyConservationTest", "LateConservationTest", "ReverseHourglassTest")))
            stop("Please choose a 'TestStatistic' that is supported by this function. E.g. TestStatistic = 'FlatLineTest', TestStatistic = 'ReductiveHourglassTest', TestStatistic = 'EarlyConservationTest', TestStatistic = 'ReverseHourglassTest'.", call. = FALSE)
            
        cat("Plot signature: '",measure, "' and test statistic: '",TestStatistic,"' running ", permutations, " permutations." )
        cat("\n")

        Stage <- NULL        
        # store transcriptome index in tibble
        if (measure == "TAI") {
                TI <-
                        tibble::tibble(Stage = names(TAI(ExpressionSet)),
                                       TI = TAI(ExpressionSet))
        }
        
        if (measure == "TDI") {
                TI <-
                        tibble::tibble(Stage = names(TDI(ExpressionSet)),
                                       TI = TDI(ExpressionSet))
        }
        
        if (measure == "TPI") {
                TI <-
                        tibble::tibble(Stage = names(TPI(ExpressionSet)),
                                       TI = TPI(ExpressionSet))
        }
        
        
        
        if (TestStatistic == "FlatLineTest") {
                if (is.null(custom.perm.matrix)) {
                        resList <- FlatLineTest(ExpressionSet = ExpressionSet,
                                                permutations  = permutations)
                }
                
                else if (!is.null(custom.perm.matrix)) {
                        resList <- FlatLineTest(ExpressionSet      = ExpressionSet,
                                                custom.perm.matrix = custom.perm.matrix)
                        
                }
        }
        
        if (TestStatistic == "ReductiveHourglassTest") {
                if (lillie.test) {
                        if (is.null(custom.perm.matrix)) {
                                resList <- ReductiveHourglassTest(
                                        ExpressionSet = ExpressionSet,
                                        modules       = modules,
                                        permutations  = permutations,
                                        lillie.test   = TRUE
                                )
                        }
                        
                        else if (!is.null(custom.perm.matrix)) {
                                resList <-
                                        ReductiveHourglassTest(
                                                ExpressionSet      = ExpressionSet,
                                                modules            = modules,
                                                lillie.test        = TRUE,
                                                custom.perm.matrix = custom.perm.matrix
                                        )
                        }
                }
                
                if (!lillie.test) {
                        if (is.null(custom.perm.matrix)) {
                                resList <- ReductiveHourglassTest(
                                        ExpressionSet = ExpressionSet,
                                        modules       = modules,
                                        permutations  = permutations,
                                        lillie.test   = FALSE
                                )
                        }
                        
                        else if (!is.null(custom.perm.matrix)) {
                                resList <-
                                        ReductiveHourglassTest(
                                                ExpressionSet      = ExpressionSet,
                                                modules            = modules,
                                                lillie.test        = FALSE,
                                                custom.perm.matrix = custom.perm.matrix
                                        )
                        }
                }
        }
        
        if (TestStatistic == "ReverseHourglassTest") {
                if (lillie.test) {
                        if (is.null(custom.perm.matrix)) {
                                resList <- ReverseHourglassTest(
                                        ExpressionSet = ExpressionSet,
                                        modules       = modules,
                                        permutations  = permutations,
                                        lillie.test   = TRUE
                                )
                        }
                        
                        else if (!is.null(custom.perm.matrix)) {
                                resList <-
                                        ReverseHourglassTest(
                                                ExpressionSet      = ExpressionSet,
                                                modules            = modules,
                                                lillie.test        = TRUE,
                                                custom.perm.matrix = custom.perm.matrix
                                        )
                        }
                }
                
                if (!lillie.test) {
                        if (is.null(custom.perm.matrix)) {
                                resList <- ReverseHourglassTest(
                                        ExpressionSet = ExpressionSet,
                                        modules       = modules,
                                        permutations  = permutations,
                                        lillie.test   = FALSE
                                )
                        }
                        
                        else if (!is.null(custom.perm.matrix)) {
                                resList <-
                                        ReverseHourglassTest(
                                                ExpressionSet      = ExpressionSet,
                                                modules            = modules,
                                                lillie.test        = FALSE,
                                                custom.perm.matrix = custom.perm.matrix
                                        )
                        }
                }
        }
        
        if (TestStatistic == "EarlyConservationTest") {
                if (lillie.test) {
                        if (is.null(custom.perm.matrix)) {
                                resList <- EarlyConservationTest(
                                        ExpressionSet = ExpressionSet,
                                        modules       = modules,
                                        permutations  = permutations,
                                        lillie.test   = TRUE
                                )
                        }
                        
                        else if (!is.null(custom.perm.matrix)) {
                                resList <-
                                        EarlyConservationTest(
                                                ExpressionSet      = ExpressionSet,
                                                modules            = modules,
                                                lillie.test        = TRUE,
                                                custom.perm.matrix = custom.perm.matrix
                                        )
                        }
                }
                
                if (!lillie.test) {
                        if (is.null(custom.perm.matrix)) {
                                resList <- EarlyConservationTest(
                                        ExpressionSet = ExpressionSet,
                                        modules       = modules,
                                        permutations  = permutations,
                                        lillie.test   = FALSE
                                )
                        }
                        
                        else if (!is.null(custom.perm.matrix)) {
                                resList <-
                                        EarlyConservationTest(
                                                ExpressionSet      = ExpressionSet,
                                                modules            = modules,
                                                lillie.test        = FALSE,
                                                custom.perm.matrix = custom.perm.matrix
                                        )
                        }
                        
                }
        }
        
        
        if (TestStatistic == "LateConservationTest") {
                     if (lillie.test) {
                        if (is.null(custom.perm.matrix)) {
                              resList <- LateConservationTest(
                                      ExpressionSet = ExpressionSet,
                                      modules       = modules,
                                      permutations  = permutations,
                                      lillie.test   = TRUE
                          )
                     }
            
                    else if (!is.null(custom.perm.matrix)) {
                              resList <-
                                      LateConservationTest(
                                        ExpressionSet      = ExpressionSet,
                                        modules            = modules,
                                        lillie.test        = TRUE,
                                        custom.perm.matrix = custom.perm.matrix
                                      )
                    }
          }
          
          if (!lillie.test) {
                    if (is.null(custom.perm.matrix)) {
                            resList <- LateConservationTest(
                                      ExpressionSet = ExpressionSet,
                                      modules       = modules,
                                      permutations  = permutations,
                                      lillie.test   = FALSE
                            )
                    }
            
                    else if (!is.null(custom.perm.matrix)) {
                            resList <-
                                    LateConservationTest(
                                      ExpressionSet      = ExpressionSet,
                                      modules            = modules,
                                      lillie.test        = FALSE,
                                      custom.perm.matrix = custom.perm.matrix
                                    )
                    }
            
          }
        }
        
        # get p-value and standard deviation values from the test statistic
        pval <- resList$p.value
        pval <- format(pval,digits = 3)
        std_dev <- resList$std.dev
        
        TI.ggplot <- ggplot2::ggplot(TI, ggplot2::aes(
                x = factor(Stage, levels = unique(Stage)),
                y = TI,
                group = 1
        ))  + ggplot2::geom_ribbon(ggplot2::aes(
                ymin = TI - std_dev,
                ymax = TI + std_dev
        ), alpha = alpha) +
                ggplot2::geom_line(lwd = lwd) +
                ggplot2::theme_classic() +
                ggplot2::labs(x = xlab, y = ylab, title = main) +
                ggplot2::theme(
                        title            = ggplot2::element_text(size = 12),
                        aspect.ratio     = 0.55,
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
                ggplot2::scale_y_continuous(breaks = scales::pretty_breaks(n = y.ticks))
        
        if ((TestStatistic == "FlatLineTest") && p.value){
                TI.ggplot <-
                        TI.ggplot + ggplot2::labs(
                                subtitle = paste0("p_flt = ", pval),
                                size = 6
                        )
                cat("\n")
                cat("Significance status of signature: ", ifelse(as.numeric(pval) <= 0.05, "significant.","not significant (= no evolutionary signature in the transcriptome)."))
        }
                
        if (TestStatistic == "ReductiveHourglassTest"){
                
                if (p.value) {
                        TI.ggplot <-
                                TI.ggplot + ggplot2::labs(
                                  subtitle =  paste0("p_rht = ", pval),
                                  size = 6
                                )  
                }
                
                if (shaded.area) {
                        TI.ggplot <-
                                TI.ggplot + ggplot2::geom_rect(data = TI,ggplot2::aes(
                                        xmin = modules[[2]][1],
                                        xmax = modules[[2]][length(modules[[2]])],
                                        ymin = min(TI) - (min(TI) / 50),
                                        ymax = Inf), fill = "#4d004b", alpha = alpha)  
                }
                
                stage.names <- names(ExpressionSet)[3:ncol(ExpressionSet)]
                cat("Modules: \n early = {",paste0(stage.names[modules[[1]]], " "),"}","\n","mid = {",paste0(stage.names[modules[[2]]], " "),"}","\n","late = {",paste0(stage.names[modules[[3]]], " "),"}")
                cat("\n")
                cat("Significance status of signature: ", ifelse(as.numeric(pval) <= 0.05, "significant.","not significant (= no evolutionary signature in the transcriptome)."))
        }   
        
        if (TestStatistic == "ReverseHourglassTest"){
                
                if (p.value) {
                        TI.ggplot <-
                                TI.ggplot + ggplot2::labs(
                                  subtitle = paste0("p_reverse_hourglass = ", pval),
                                  size = 6
                                )  
                }
                
                if (shaded.area) {
                        TI.ggplot <-
                                TI.ggplot + ggplot2::geom_rect(data = TI,ggplot2::aes(
                                        xmin = modules[[2]][1],
                                        xmax = modules[[2]][length(modules[[2]])],
                                        ymin = min(TI) - (min(TI) / 50),
                                        ymax = Inf), fill = "#4d004b", alpha = alpha)  
                }
                
                stage.names <- names(ExpressionSet)[3:ncol(ExpressionSet)]
                cat("Modules: \n early = {",paste0(stage.names[modules[[1]]], " "),"}","\n","mid = {",paste0(stage.names[modules[[2]]], " "),"}","\n","late = {",paste0(stage.names[modules[[3]]], " "),"}")
                cat("\n")
                cat("Significance status of signature: ", ifelse(as.numeric(pval) <= 0.05, "significant.","not significant (= no evolutionary signature in the transcriptome)."))
        }  
        
        if (TestStatistic == "EarlyConservationTest"){
                
                if (p.value) {
                        TI.ggplot <-
                                TI.ggplot + ggplot2::labs(
                                  subtitle = paste0("p_ect = ", pval),
                                  size = 6
                                )  
                }
                
                if (shaded.area) {
                        TI.ggplot <-
                                TI.ggplot + ggplot2::geom_rect(data = TI,ggplot2::aes(
                                        xmin = modules[[2]][1],
                                        xmax = modules[[2]][length(modules[[2]])],
                                        ymin = min(TI) - (min(TI) / 50),
                                        ymax = Inf), fill = "#4d004b", alpha = alpha * 0.5)  
                }
                
                stage.names <- names(ExpressionSet)[3:ncol(ExpressionSet)]
                cat("Modules: \n early = {",paste0(stage.names[modules[[1]]], " "),"}","\n","mid = {",paste0(stage.names[modules[[2]]], " "),"}","\n","late = {",paste0(stage.names[modules[[3]]], " "),"}")
                cat("\n")
                cat("Significance status of signature: ", ifelse(as.numeric(pval) <= 0.05, "significant.","not significant (= no evolutionary signature in the transcriptome)."))
        }
        
        if (TestStatistic == "LateConservationTest"){
          
                if (p.value) {
                        TI.ggplot <-
                          TI.ggplot + ggplot2::labs(
                            subtitle = paste0("p_lct = ", pval),
                            size = 6
                          )  
                }
                
                if (shaded.area) {
                        TI.ggplot <-
                          TI.ggplot + ggplot2::geom_rect(data = TI,ggplot2::aes(
                            xmin = modules[[2]][1],
                            xmax = modules[[2]][length(modules[[2]])],
                            ymin = min(TI) - (min(TI) / 50),
                            ymax = Inf), fill = "#4d004b", alpha = alpha * 0.5)  
                }
                
                stage.names <- names(ExpressionSet)[3:ncol(ExpressionSet)]
                cat("Modules: \n early = {",paste0(stage.names[modules[[1]]], " "),"}","\n","mid = {",paste0(stage.names[modules[[2]]], " "),"}","\n","late = {",paste0(stage.names[modules[[3]]], " "),"}")
                cat("\n")
                cat("\n")
                cat("Significance status of signature: ", ifelse(as.numeric(pval) <= 0.05, "significant.","not significant (= no evolutionary signature in the transcriptome)."))
        }
        
        TI.ggplot <- TI.ggplot + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 1,hjust = 1))
        
        if (TestStatistic == "FlatLineTest")
                message("\n")
                message("-> Now run 'FlatLineTest(..., permutations  = ", permutations, ", plotHistogram = TRUE)' to analyse the permutation test performance.")
        return (TI.ggplot)
}
