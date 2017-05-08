




## Define function ####
`rrarefy.record` <-
  function (counts.file, site.name, zones=NULL, n.random=1000)
  {
    
    ## Load packages
    require(vegan) # to use function rrarefy()
    require(scales) # for transparency of layers
    
    # Define and create output directory
    output.dir <- "./Out/"
    dir.create(output.dir)
    
    # Define input data and input parameter
    x <- counts.file[ ,3:ncol(counts.file)]
    scale.depth <- counts.file[ ,1]
    scale.age <- counts.file[ ,2]
    x.depth = c(max(scale.depth), min(scale.depth))
    x.age = c(max(scale.age), min(scale.age))
    
    # Check, if no zones defined in function, then uses top and bottom age limits
    if(is.null(zones)) {
      zones <- x.age
    } else { }
    
    
    
    ## Calculate rarefied richness and diversity (N0, N1, N2, etc.) ####
    
    # Make space in empty data.frame 'y' where output will be stored
    raremax <- min(rowSums(x))
    n.row <- nrow(x)
    y <- matrix(data=0, nrow=n.row, ncol=12)
    y <- as.data.frame(y)
    colnames(y) <- c("N0.lci","N0","N0.uci","N1.lci","N1","N1.uci","N2.lci","N2","N2.uci",
                     "N2.N0.lci","N2.N0","N2.N0.uci")
    
    # Make space in matrices, one for each variable to be calculated
    N0.r <- matrix(data=NA, nrow=n.row, ncol=n.random)
    N1.r <- matrix(data=NA, nrow=n.row, ncol=n.random)
    N2.r <- matrix(data=NA, nrow=n.row, ncol=n.random)
    N2.N0.r  <- matrix(data=NA, nrow=n.row, ncol=n.random)
    N2.N1.r <- matrix(data=NA, nrow=n.row, ncol=n.random)
    
    ## Take random sample from each Sample 1000 times,
    ## calculate rarefied richness (N0) and rarefied diversity (N1 and N2), and
    ## populate matrices prepared above
    for (i in 1:n.random) {
      rand <- rrarefy(x, sample=raremax)
      N0.i <- specnumber(rand)
      N0.r[, i] <- N0.i
      
      N1.i <- exp(diversity(rand, index="shannon"))
      N1.r[, i] <- N1.i
      
      N2.i <- diversity(rand, index="invsimpson")
      N2.r[, i] <- N2.i
      
      N2.N0.r[, i] <- N2.i / N0.i
      N2.N1.r[, i] <- N2.i / N1.i
    }
    
    rm(i, rand, N0.i, N1.i, N2.i)
    
    
    # Extract means and sd from rarefied richness and diversity
    # and populate data.frame 'y' that was prepared above
    y$N0 <- apply(X=N0.r, MARGIN=1, FUN=mean)
    st.dev <- apply(X=N0.r, MARGIN=1, FUN=sd)
    y$N0.uci <- qnorm(0.975, mean=y$N0, sd=st.dev)
    y$N0.lci <- qnorm(0.975, mean=y$N0, sd=st.dev, lower.tail=FALSE)
    
    y$N1 <- apply(X=N1.r, MARGIN=1, FUN=mean)
    st.dev <- apply(X=N1.r, MARGIN=1, FUN=sd)
    y$N1.uci <- qnorm(0.975, mean=y$N1, sd=st.dev)
    y$N1.lci <- qnorm(0.975, mean=y$N1, sd=st.dev, lower.tail=FALSE)
    
    y$N2 <- apply(X=N2.r, MARGIN=1, FUN=mean)
    st.dev <- apply(X=N2.r, MARGIN=1, FUN=sd)
    y$N2.uci <- qnorm(0.975, mean=y$N2, sd=st.dev)
    y$N2.lci <- qnorm(0.975, mean=y$N2, sd=st.dev, lower.tail=FALSE)
    
    
    y$N2.N0 <- apply(X=N2.N0.r, MARGIN=1, FUN=mean)
    st.dev <- apply(X=N2.N0.r, MARGIN=1, FUN=sd)
    y$N2.N0.uci <- qnorm(0.975, mean=y$N2.N0, sd=st.dev)
    y$N2.N0.lci <- qnorm(0.975, mean=y$N2.N0, sd=st.dev, lower.tail=FALSE)
    
    y$N2.N1 <- apply(X=N2.N1.r, MARGIN=1, FUN=mean)
    
    rm(N0.r, N1.r, N2.r, st.dev, n.random)
    
    ## Add depth and age scales, and save data.frame 'y' as .csv file in output directory
    y$depth <- scale.depth
    y$age <- scale.age
    y <- y[ ,c(14,15,1:13)]
    
    write.csv(y, file=paste0(output.dir, "richness_diversity.csv"), row.names=F)
    
    
    ## Make plots ####
    pdf(file.path(output.dir, "Figure 1.pdf"))
    par(mfrow=c(2,1), mar=c(0.2,3.5,1,0), oma=c(5,0.5,0.5,2))
    plot(y$age, y$N0, xlim=x.age, ylim=c(5,35), type="l", ylab="", axes=F, col="red", lwd=2, main=site.name)
    polygon(x=c(y$age, rev(y$age)), y=c(y$N0.uci, rev(y$N0.lci)), col=alpha("red", 0.3), border=F)
    lines(y$age, y$N1, xlim=x.age, type="l", ylab="N1", col="blue", lwd=2)
    polygon(x=c(y$age, rev(y$age)), y=c(y$N1.uci, rev(y$N1.lci)), col=alpha("blue", 0.3), border=F)
    lines(y$age, y$N2, xlim=x.age, type="l", ylab="N2", col="green", lwd=2)
    polygon(x=c(y$age, rev(y$age)), y=c(y$N2.uci, rev(y$N2.lci)), col=alpha("green", 0.3), border=F)
    abline(v=zones, col="grey")
    axis(1, labels=F)
    axis(2, las=2)
    mtext("N", side=2, line=3)
    legend(x="topleft", inset=c(0,0), legend=c("N0","N1","N2"), fill=c("red","blue","green"),
           border="NA", bty="n", horiz=T, cex=0.9)
    
    plot(y$age, y$N2.N0, xlim=x.age, type="l", ylab="", axes=F, col="red", lwd=2)
    polygon(x=c(y$age, rev(y$age)), y=c(y$N2.N0.uci, rev(y$N2.N0.lci)), col=alpha("red", 0.3), border=F)
    axis(1, labels=F)
    axis(2, las=2)
    mtext("N2/N0", side=2, line=3)
    abline(v=zones, col="grey")
    axis(1, labels=T)
    mtext("Age (cal yrs BP)", side=1, line=2)
    dev.off()
    
  }
