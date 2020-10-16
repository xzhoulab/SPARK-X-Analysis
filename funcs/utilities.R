## QQ Plots
##==============
library(tidyr)
library(ggplot2)
library(reshape2)
library(scales) 
qplot_gg <- function(pvalues,legend.position="top",ylab=NULL,xlab=NULL,ax.txt.size=20,ax.title.size=20,len.title.size=1.5,len.txt.size=1,len.title=NULL,col.alpha=1,pt.size=2,
    col.base=NULL,cl=1,single_col=NULL,keeporder=T,factor_level=NULL,self_label=NULL){

    if(is.null(col.base)){
        if(cl==1){
            # col.base <- c("hotpink","mediumorchid2","#F4A582","#B8E186","#8dd3c7","wheat")
            col.base <- c("hotpink","mediumorchid2","lightskyblue")
        }else{
            if(is.list(pvalues)){
                col.base <- hue_pal()(length(pvalues))
            }else{
                col.base <- hue_pal()(10)
            } 
        }
    }

    if(is.null(ylab)){ylab=expression(paste("Observed ",-log[10],"(",italic(p),"-value)","\n"))}
    if(is.null(xlab)){xlab=expression(paste("Expected ",-log[10],"(",italic(p),"-value)"))}
    if(is.list(pvalues)){
        nn          <- length(pvalues) 
        exp.vec     <- c()
        for(iset in 1:nn){
            n           <- length(pvalues[[iset]])
            exp.x       <- -log10((rank(pvalues[[iset]], ties.method="first")-.5)/n)
            exp.vec     <- c(exp.vec,exp.x)
        }

        pd              <- melt(pvalues)
        pd$pve          <- exp.vec
        pd$pvo          <- -log10(pd$value)

        if(keeporder){
            if(is.null(factor_level)){
                factor_level <- unique(pd$L1)
            }
            pd$L1       <- factor(pd$L1,levels=factor_level)
        }

        n = length(unique(exp.vec))
        df1 <- data.frame(
                expected = -log10(ppoints(n)),
                clower   = -log10(qbeta(p = (1 - 0.95) / 2, shape1 = 1:n, shape2 = n:1)),
                cupper   = -log10(qbeta(p = (1 + 0.95) / 2, shape1 = 1:n, shape2 = n:1))
                )
        p1 <- ggplot(df1) +geom_ribbon(aes(x=expected, ymax=cupper, ymin=clower), fill="darkgrey", alpha=.5)+
                geom_abline(intercept = 0, slope = 1, size = 1.5,col="white")
        p2 <- p1 + geom_point(data=pd, aes(pve, pvo,color = L1), size = pt.size,alpha = col.alpha) +
                scale_color_manual(values=col.base,labels = self_label)+
                scale_x_continuous(xlab) +
                scale_y_continuous(ylab) +
                labs(color=len.title)+
                theme_bw()+
                theme(axis.text.x = element_text(size=ax.txt.size),axis.text.y=element_text(size=ax.txt.size),axis.title.y = element_text(size=ax.title.size), axis.title.x = element_text(size=ax.title.size),
                        legend.text= element_text(size=rel(len.txt.size)),legend.title = element_text(size=rel(len.title.size)),
                        panel.border = element_blank(), panel.grid.major = element_blank(),
                        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),legend.position=legend.position,legend.text.align = 0)
    
    }else{
        n   <- length(pvalues)
        df1 <- data.frame(
            expected = -log10(ppoints(n)),
            observed = -log10(sort(pvalues)),
            clower   = -log10(qbeta(p = (1 - 0.95) / 2, shape1 = 1:n, shape2 = n:1)),
            cupper   = -log10(qbeta(p = (1 + 0.95) / 2, shape1 = 1:n, shape2 = n:1))
        )

        if(is.null(single_col)){single_col = "lightpink"}

        p1  <- ggplot(df1) +geom_ribbon(aes(x=expected, ymax=cupper, ymin=clower), fill="darkgrey", alpha=.5)+
                geom_abline(intercept = 0, slope = 1, size = 1.5,col="white")+
                geom_point(aes(expected, observed), size = 3,col=single_col,alpha = col.alpha) 

        p2  <- p1 + scale_x_continuous(xlab) +
                scale_y_continuous(ylab) +
                theme_bw()+
                theme(axis.text.x = element_text(size=ax.txt.size),axis.text.y=element_text(size=ax.txt.size),
                        axis.title.y = element_text(size=ax.title.size), axis.title.x = element_text(size=ax.title.size),
                        legend.text= element_text(size=rel(len.txt.size)),legend.title = element_text(size=rel(len.title.size)),
                        panel.border = element_blank(), panel.grid.major = element_blank(),
                        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),legend.position=legend.position)
    }
    return(p2)
}




## Calculate Empirical FDR
##============================

fdr.cal <- function(altpval,nullpval,rn=NULL){ 
    if(!is.null(rn)){
            order.rn <- rn[order(altpval)]
        }else{
            order.rn <- NULL
        }
    altpval.order   <- altpval[order(altpval)]
    numperm         <- ncol(nullpval)
    nullpval.order  <- sapply(1:numperm,function(x){nullpval[,x][order(nullpval[,x])]})
    out             <- c()
    for(i in 1:nrow(nullpval)){
        s = sum(nullpval.order[1:i,]<=altpval.order[i])
        e   <- s/(i*numperm)
        out <- c(out,e)
    }
    if(is.null(order.rn)){
        return(out)
    }else{
        res <- cbind.data.frame(rn=order.rn,fdr=out)
        return(res)
    }
    
}

