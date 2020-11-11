
##############################################################################################
#  Functions to extract estimates from the summary table and from the posterior sim matrix

# Extract effect esttimates from the summary table
extract.effect <- function(fit, effect.name, index.type=NA, source="Bayesian"){
  # index.type must be one of "year", "site", or "year.site" or "year.subregion"
  #browser()
  select <- grepl(effect.name, rownames(fit$results$BUGSoutput$summary))
  effect <- as.data.frame(fit$results$BUGSoutput$summary[select,c("mean", "sd", "2.5%","97.5%","Rhat", "n.eff"),drop=FALSE])
  names(effect) <- make.names(names(effect))
  if("year" == index.type){
     effect$year.index <- as.numeric(substr(rownames(effect),
                                          1+regexpr("[", rownames(effect), fixed=TRUE),
                                         -1+regexpr("]", rownames(effect), fixed=TRUE)))
     effect <- merge(effect, fit$year.index)
  }
  if("site" == index.type){
     effect$site.index <- as.numeric(substr(rownames(effect),
                                          1+regexpr("[", rownames(effect), fixed=TRUE),
                                         -1+regexpr("]", rownames(effect), fixed=TRUE)))
     effect <- merge(effect, fit$site.index)
  }
  if("year.site"== index.type){
    effect$year.index <- as.numeric(substr(rownames(effect),
                                          1+regexpr("[", rownames(effect), fixed=TRUE),
                                         -1+regexpr(",", rownames(effect), fixed=TRUE)))
    effect$site.index <- as.numeric(substr(rownames(effect),
                                          1+regexpr(",", rownames(effect), fixed=TRUE),
                                         -1+regexpr("]", rownames(effect), fixed=TRUE)))
    effect <- merge(effect, fit$site.index)
    effect <- merge(effect, fit$year.index)
  }
  if("year.subregion"==index.type){
    effect$year.index <- as.numeric(substr(rownames(effect),
                                          1+regexpr("[", rownames(effect), fixed=TRUE),
                                         -1+regexpr(",", rownames(effect), fixed=TRUE)))
    effect$subregion.index <- as.numeric(substr(rownames(effect),
                                          1+regexpr(",", rownames(effect), fixed=TRUE),
                                         -1+regexpr("]", rownames(effect), fixed=TRUE)))
    effect <- merge(effect, fit$year.index)
    effect <- merge(effect, fit$subregion.index)
  }
  effect$Source <- source
  effect
}



extract.posterior <- function(fit, effect.name, index.type=" ", source="Bayesian"){
  # extract the posterior for an effect in the long-format
  # index.type must be one of "year", "site", or "year.site"
  #browser()
  select <- grepl(effect.name, colnames(fit$results$BUGSoutput$sims.matrix) )
  post   <- as.data.frame(fit$results$BUGSoutput$sims.matrix[,select,drop=FALSE])

  post$sim <- 1:nrow(post)
  post.long <- reshape2::melt(post,
                              id.vars="sim",
                              value.name="value")
  #browser()
  # If there is one index, it is could be a site or a year
  if("year" == index.type){
     post.long$year.index <- as.numeric(substr(post.long$variable,
                                          1+regexpr("[", post.long$variable, fixed=TRUE),
                                         -1+regexpr("]", post.long$variable, fixed=TRUE)))
     post.long <- merge(post.long, fit$year.index)
  }
  if("site" == index.type){
     post.long$site.index <- as.numeric(substr(post.long$variable,
                                          1+regexpr("[", post.long$variable, fixed=TRUE),
                                         -1+regexpr("]", post.long$variable, fixed=TRUE)))
     post.long <- merge(post.long, fit$site.index)
  }
  if("year.site"== index.type){
    post.long$year.index <- as.numeric(substr(post.long$variable,
                                          1+regexpr("[", post.long$variable, fixed=TRUE),
                                         -1+regexpr(",", post.long$variable, fixed=TRUE)))
    post.long$site.index <- as.numeric(substr(post.long$variable,
                                          1+regexpr(",", post.long$variable, fixed=TRUE),
                                         -1+regexpr("]", post.long$variable, fixed=TRUE)))
    post.long <- merge(post.long, fit$site.index)
    post.long <- merge(post.long, fit$year.index)
  }
  post.long$Source <- source
  post.long$variable <- as.character(post.long$variable)
  post.long
}


# Extract the tukey-freeman statistics
extract.TF.post <- function(fit){
  require(reshape2)
  select <- grepl("^TF.", colnames(fit$results$BUGSoutput$sims.matrix) )
  TF.post <- as.data.frame(fit$results$BUGSoutput$sims.matrix[,select])
  TF.post
}


postPriorOverlap <- function (paramSampleVec, prior, ..., yaxt = "n", ylab = "", 
    xlab = "Parameter", main = "", cex.lab = 1.5, cex = 1.4, 
    xlim = range(paramSampleVec), breaks = NULL, plot=FALSE) 
{
  # This is a copy of the postPriorOverlap function from the BEST package
  # except that we stopped the display of the histograms on the plot window
  require(HDInterval)
    oldpar <- par(xpd = NA)
    on.exit(par(oldpar))
    if (is.null(breaks)) {
        nbreaks <- ceiling(diff(range(paramSampleVec))/as.numeric(diff(hdi(paramSampleVec))/18))
        breaks <- seq(from = min(paramSampleVec), to = max(paramSampleVec), 
            length.out = nbreaks)
    }
    histinfo <- hist(paramSampleVec, xlab = xlab, yaxt = yaxt, 
        ylab = ylab, freq = FALSE, border = "white", col = "skyblue", 
        xlim = xlim, main = main, cex = cex, cex.lab = cex.lab, 
        breaks = breaks, plot=plot)
    if (is.numeric(prior)) {
        priorInfo <- hist(prior, breaks = c(-Inf, breaks, Inf), 
            add = TRUE, freq = FALSE, col = "yellow", border = "white", plot=plot)$density[2:length(breaks)]
    }
    else if (is.function(prior)) {
        if (class(try(prior(0.5, ...), TRUE)) == "try-error") 
            stop(paste("Incorrect arguments for the density function", 
                substitute(prior)))
        priorInfo <- prior(histinfo$mids, ...)
    }
    minHt <- pmin(priorInfo, histinfo$density)
    if(plot){
      rect(breaks[-length(breaks)], rep(0, length(breaks) - 1), 
        breaks[-1], minHt, col = "green", border = "white")
    }
    overlap <- sum(minHt * diff(histinfo$breaks))
    if (is.function(prior)){
        if(plot){
          lines(histinfo$mids, priorInfo, lwd = 2, col = "brown")
          text(mean(breaks), 0, paste0("overlap = ", round(overlap * 100), "%"), pos = 3, cex = cex)
        }
    }
    return(overlap)
}

###########################################################################
# take a data frame of estimates, 95% ci, select two sources and plot

compare.est.plot <- function(df, x.source, y.source, source.name="Source", 
                             est.name="mean", lcl.name="X2.5.", ucl.name="X97.5.",
                             year.name="year"){
    # df must have columns for mean, x2.5., X97.5
    df <- df[ df[,source.name] %in% c(x.source, y.source),]  # select the two variables
    df.wide.est <- reshape2::dcast(df, as.formula(paste(year.name,"~",source.name)), value.var=est.name)
    plot <- ggplot(data=df.wide.est, aes_string(x=x.source, y=y.source))+
       geom_point()+
       geom_abline(intercept=0, slope=1)+
       geom_smooth(method="lm", se=FALSE)+
       geom_text(aes_string(label=year.name))
    # get the confidence limits for the x-variable
    x.ci <- df[ df[,source.name] %in% x.source,]
    df.wide.xci <- merge(df.wide.est, x.ci[,c(year.name, lcl.name, ucl.name)], by=year.name)
    plot <- plot +
       geom_errorbarh(data=df.wide.xci, aes_string(xmin=lcl.name, xmax=ucl.name), alpha=0.5, width=0)
    # get the confidence limits for the y-variable
    y.ci <- df[ df[,source.name] %in% y.source,]
    df.wide.yci <- merge(df.wide.est, y.ci[,c(year.name, lcl.name, ucl.name)], by=year.name)
    plot <- plot +
       geom_errorbar(data=df.wide.yci, aes_string(ymin=lcl.name, ymax=ucl.name), alpha=0.5, width=0)
    plot <- plot +
           xlab(paste0(x.source," (95% ci)"))+
           ylab(paste0(y.source," (95% ci)"))
           #+ xlim(0,1)+ ylim(0,1)
    plot
}

#compare.est.plot(plotdata, "LSMeans","MM.p.uw")


