##
##  BioCro/R/c4photo.R by Fernando Ezequiel Miguez  Copyright (C) 2007-2008
##
##  This program is free software; you can redistribute it and/or modify
##  it under the terms of the GNU General Public License as published by
##  the Free Software Foundation; either version 2 or 3 of the License
##  (at your option).
##
##  This program is distributed in the hope that it will be useful,
##  but WITHOUT ANY WARRANTY; without even the implied warranty of
##  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##  GNU General Public License for more details.
##
##  A copy of the GNU General Public License is available at
##  http://www.r-project.org/Licenses/
##
##

c4photo <- function(Qp,Tl,RH,vmax=39,alpha=0.04,kparm=0.7,theta=0.83,
                    beta=0.93,Rd=0.8,Catm=380,b0=0.08,b1=3,
                    stress=1,ws=c("gs","vmax"))
{
    if((max(RH) > 1) || (min(RH) < 0))
        stop("RH should be between 0 and 1")
    if(any(Catm < 150))
        warning("Stomatal conductance is not reliable for values of Catm lower than 150\n")
    if(any(Catm < 15))
        warning("Assimilation is not reliable for low (<15) Catm values")
    ws <- match.arg(ws)
    if(ws == "gs") ws <- 1
    else ws <- 0

    if(length(Catm) == 1){
      Catm <- rep(Catm,length(Qp))
    }else{
      if(length(Catm) != length(Qp))
        stop("length of Catm should be either 1 or equal to length of Qp")
    }
    
    res <- .Call(c4photo_sym,as.double(Qp),
                 as.double(Tl),as.double(RH),
                 as.double(vmax),as.double(alpha),
                 as.double(kparm),as.double(theta),
                 as.double(beta),
                 as.double(Rd),as.double(Catm),
                 as.double(b0),as.double(b1),as.double(stress),as.integer(ws))
    res
}


MCMCc4photo <- function(data, niter = 20000, op.level=1, ivmax = 39,
                        ialpha = 0.04, ikparm = 0.7, itheta=0.83,
                        ibeta=0.93, iRd = 0.8, Catm = 380,
                        b0 = 0.08, b1 = 3, stress=1, ws=c("gs","vmax"), scale = 1,
                        sds=c(1,0.005,0.5),prior=c(39,10,0.04,0.02,3,1.5)){

    if(ncol(data) != 4)
        stop("ncol data should be 4")

    if(op.level == 2 && length(prior) != 6){
            stop("length of prior should be 6")
    }
    
    if(niter < 2)
        stop("niter should be at least 2")
    
    assim <- data[,1]
    qp <- data[,2]
    temp <- data[,3]
    rh <- data[,4]

    ws1 <- ws
    ws <- match.arg(ws)
    if(ws == "gs") ws <- 1
    else ws <- 0

    if(op.level != 1 && op.level != 2) stop("op.level should be either 1 or 2")
    
    res <- .Call(McMCc4photo, as.double(assim), as.double(qp),
                 as.double(temp), as.double(rh), as.integer(niter),
                 as.double(ivmax), as.double(ialpha), as.double(ikparm),
                 as.double(itheta), as.double(ibeta),
                 as.double(iRd), as.double(Catm), as.double(b0), as.double(b1),
                 as.double(stress), as.double(scale), as.double(sds),
                 as.integer(ws), as.double(prior), as.integer(op.level))
    res$resuMC <- t(res$resuMC)
    res$op.level <- op.level
    res$niter <- niter
    colnames(res$resuMC) <- c("Vcmax","Alpha","Rd","RSS")
    res$prior <- prior
    res$obs <- data
    res$xpar <- list(kparm = ikparm, theta = itheta, beta = ibeta,
                     Catm = Catm, b0 = b0, b1 = b1, stress = stress, ws=ws1)
    structure(res, class = "MCMCc4photo")
}




## Function for printing the MCMCc4photo objects

print.MCMCc4photo <- function(x,burnin=1,level=0.95,digits=1,...){

    op.level <- x$op.level
    ul <- 1 - (1-level)/2
    ll <- (1 - level)/2
    if(op.level == 1){
        xMat <- x$resuMC[burnin:x$niter,1:2]
        colnames(xMat) <- c("Vmax","alpha")
    }else{
        xMat <- x$resuMC[burnin:x$niter,1:3]
        colnames(xMat) <- c("Vmax","alpha","Rd")
    }
    cat("\n Markov chain Monte Carlo for the Collatz C4 photosynthesis model")
    
    cat("\n Summary of the chain")
    cat("\n Moves:",x$accept,"Prop:",x$accept/x$niter,"\n")
    if(op.level == 1){
        cat("\n Summaries for vmax and alpha:\n")
    }else{
        cat("\n Summaries for vmax, alpha and rd:\n")
    }
    sum1 <- summary(x$resuMC[burnin:x$niter,1])
    sum2 <- summary(x$resuMC[burnin:x$niter,2])
    if(op.level == 2) sum3 <- summary(x$resuMC[burnin:x$niter,3])
    nm <- names(sum1)
    if(op.level == 1) mat <- matrix(rbind(sum1,sum2),nrow=2,ncol=6)
    if(op.level == 2) mat <- matrix(rbind(sum1,sum2,sum3),nrow=3,ncol=6)
    colnames(mat) <- nm
    if(op.level == 1) rownames(mat) <- c("vmax","alpha")
    if(op.level == 2) rownames(mat) <- c("vmax","alpha","rd")
    print(mat,...)
    if(op.level == 1) cat("\n",level*100,"% Quantile Intervals for vmax and alpha:\n")
    if(op.level == 2) cat("\n",level*100,"% Quantile Intervals for vmax, alpha and rd:\n")
    qua1 <- quantile(x$resuMC[burnin:x$niter,1],c(ll,ul))
    qua2 <- quantile(x$resuMC[burnin:x$niter,2],c(ll,ul))
    if(op.level == 2) qua3 <- quantile(x$resuMC[burnin:x$niter,3],c(ll,ul))
    if(op.level == 1) mat2 <- rbind(qua1,qua2)
    if(op.level == 2) mat2 <- rbind(qua1,qua2,qua3)
    if(op.level == 1) rownames(mat2) <- c("vmax","alpha")
    if(op.level == 2) rownames(mat2) <- c("vmax","alpha","rd")
    colnames(mat2) <- c(ll,ul)
    print(mat2,...)
    cat("\n correlation matrix:\n")
    print(cor(xMat),...)
    cat("\n RSS range:",range(x$resuMC[burnin:x$niter,4]),"\n")
    invisible(x)
    
}


plot.MCMCc4photo <- function(x,x2=NULL,x3=NULL,
                             plot.kind=c("trace","density","OandF"),type=c("l","p"),
                             burnin=1,cols=c("blue","green","purple"),prior=FALSE,pcol="black",...){

    plot.kind <- match.arg(plot.kind)
    type <- match.arg(type)
    op.level <- x$op.level
    ## This first code is to plot the first object only
    ## Ploting the trace
    if(missing(x2) && missing(x3)){
        if(plot.kind == "trace"){
            plot1 <-  xyplot(x$resuMC[burnin:x$niter,1] ~ burnin:x$niter ,
                             xlab = "Iterations", type = type, col=cols[1],
                             ylab = expression(paste("Vmax (",mu,mol," ",m^-2," ",s^-1,")")),
                             ...)
            plot2 <-  xyplot(x$resuMC[burnin:x$niter,2] ~ burnin:x$niter ,
                             xlab = "Iterations", type = type, col=cols[1],
                             ylab = expression(paste("alpha (",mol," ",m^-1,")")),
                             ...)
            if(op.level == 2){
                plot3 <-  xyplot(x$resuMC[burnin:x$niter,3] ~ burnin:x$niter ,
                                 xlab = "Iterations", type = type, col=cols[1],
                                 ylab = expression(paste("Rd (",mu,mol," ",m^-2," ",s^-1,")")),
                                 ...)
            }
            if(op.level == 1){
                print(plot1,position=c(0,0,0.5,1),more=TRUE)
                print(plot2,position=c(0.5,0,1,1))
            }else{
                print(plot1,position=c(0,0,0.3333,1),more=TRUE)
                print(plot2,position=c(0.3333,0,0.6666,1),more=TRUE)
                print(plot3,position=c(0.6666,0,1,1))
            }
        } else
        ## Ploting the density
        if(plot.kind == "density"){
            if(prior == FALSE){
                plot1 <-  densityplot(~x$resuMC[burnin:x$niter,1],xlab="Vmax",
                                      col=cols[1],
                                      plot.points=FALSE,...)
                plot2 <-  densityplot(~x$resuMC[burnin:x$niter,2],xlab="alpha",
                                      col=cols[1],
                                      plot.points=FALSE,...)
                if(op.level == 2){
                    plot3 <-  densityplot(~x$resuMC[burnin:x$niter,3],xlab="rd",
                                          col=cols[1],
                                          plot.points=FALSE,...)
                }
            }else{
                plot1 <-  densityplot(~x$resuMC[burnin:x$niter,1],xlab="Vmax",
                                      col=cols[1],
                                      plot.points=FALSE,
                                      panel = function(xi,...){
                                          panel.densityplot(xi,...)
                                          panel.mathdensity(dmath=dnorm, args=list(mean = x$prior[1], sd = x$prior[2]), col=pcol)
                                      },...)
                plot2 <-  densityplot(~x$resuMC[burnin:x$niter,2],xlab="alpha",
                                      col=cols[1],
                                      plot.points=FALSE,
                                      panel = function(xi,...){
                                          panel.densityplot(xi,...)
                                          panel.mathdensity(dmath=dnorm, args=list(mean = x$prior[3], sd = x$prior[4]), col=pcol)
                                      },...)
                if(op.level == 2){
                    plot3 <-  densityplot(~x$resuMC[burnin:x$niter,3],xlab="rd",
                                          col=cols[1],
                                          plot.points=FALSE,
                                          panel = function(xi,...){
                                              panel.densityplot(xi,...)
                                              panel.mathdensity(dmath=dnorm, args=list(mean = x$prior[5], sd = x$prior[6]), col=pcol)
                                          },...)
                }
            }
            if(op.level == 1){
                print(plot1,position=c(0,0,0.5,1),more=TRUE)
                print(plot2,position=c(0.5,0,1,1))
            }else{
                print(plot1,position=c(0,0,0.3333,1),more=TRUE)
                print(plot2,position=c(0.3333,0,0.6666,1),more=TRUE)
                print(plot3,position=c(0.6666,0,1,1))
            }
        }
        if(plot.kind == "OandF"){
            if(!missing(x2)) stop("This option only works for one object")
            if(x$niter > 1e4) stop("Too slow for this number of iterations")
            obs <- x$obs
            xpar <- x$xpar
            obs.o <- obs[order(obs[,2]),]
            prds <- predict(x, obs = obs.o, burnin = burnin,
                            kparm = xpar$kparm, theta = xpar$theta,
                            beta = xpar$beta, Catm = xpar$Catm,
                            b0 = xpar$b0, b1 = xpar$b1, stress = xpar$stress,
                            ws = xpar$ws)
            plt <- xyplot(assim ~ qp, data = prds, col = "grey",
                          ylab = "CO2 uptake",
                          xlab = "Quantum flux",
                          groups = niter,
                          panel = function(x,y,...){
                              panel.xyplot(x,y,type = 'l',...)
                              panel.xyplot(x,y,type = 'a')
                              panel.xyplot(x = obs.o[,2], y = obs.o[,1], col = "red", type = "o")
                          })
            plot(plt)
        }
    } else
    ## This part of the code is to plot objects x and x2
    ## Ploting the trace
    if(missing(x3)){
        n1 <- x$niter
        n2 <- x2$niter
        maxchainLength <- max(n1,n2)
        tmpvec11 <- x$resuMC[burnin:n1,1]
        tmpvec12 <- x2$resuMC[burnin:n2,1]
        tmpvec21 <- x$resuMC[burnin:n1,2]
        tmpvec22 <- x2$resuMC[burnin:n2,2]
        if(op.level == 2){
            tmpvec31 <- x$resuMC[burnin:n1,3]
            tmpvec32 <- x2$resuMC[burnin:n2,3]
            ymin3 <- min(c(tmpvec31,tmpvec31))*0.95
            ymax3 <- max(c(tmpvec32,tmpvec32))*1.05
        }
        ymin1 <- min(c(tmpvec11,tmpvec12))*0.95
        ymax1 <- max(c(tmpvec11,tmpvec12))*1.05
        ymin2 <- min(c(tmpvec21,tmpvec22))*0.95
        ymax2 <- max(c(tmpvec21,tmpvec22))*1.05
        if(plot.kind == "trace"){
            plot1 <-  xyplot(tmpvec11 ~ burnin:n1 ,
                             xlim=c(I(burnin-0.05*maxchainLength),I(maxchainLength*1.05)),
                             ylim=c(ymin1,ymax1),
                             xlab = "Iterations", type = "l",
                             ylab = expression(paste("Vmax (",mu,mol," ",m^-2," ",s^-1,")")),
                             panel = function(x,y,...){
                                 panel.xyplot(x,y,col=cols[1],...)
                                 panel.xyplot(burnin:n2,tmpvec12,col=cols[2],...)
                             },...)
            plot2 <-  xyplot(tmpvec21 ~ burnin:n2 ,
                             xlim=c(I(burnin-0.05*maxchainLength),I(maxchainLength*1.05)),
                             ylim=c(ymin2,ymax2),
                             xlab = "Iterations", type = "l",
                             ylab = expression(paste("alpha (",mol," ",m^-1,")")),
                             panel = function(x,y,...){
                                 panel.xyplot(x,y,col=cols[1],...)
                                 panel.xyplot(burnin:n2,tmpvec22,col=cols[2],...)
                             },...)
            if(op.level == 2){
                plot3 <-  xyplot(tmpvec31 ~ burnin:n2 ,
                                 xlim=c(I(burnin-0.05*maxchainLength),I(maxchainLength*1.05)),
                                 ylim=c(ymin3,ymax3),
                                 xlab = "Iterations", type = "l",
                                 ylab = expression(paste("Rd (",mu,mol," ",m^-2," ",s^-1,")")),
                                 panel = function(x,y,...){
                                     panel.xyplot(x,y,col=cols[1],...)
                                     panel.xyplot(burnin:n2,tmpvec32,col=cols[2],...)
                                 },...)
            }
            if(op.level == 1){
                print(plot1,position=c(0,0,0.5,1),more=TRUE)
                print(plot2,position=c(0.5,0,1,1))
            }else{
                print(plot1,position=c(0,0,0.3333,1),more=TRUE)
                print(plot2,position=c(0.3333,0,0.6666,1),more=TRUE)
                print(plot3,position=c(0.6666,0,1,1))
            }
        } else
        ## ploting the density
        if(plot.kind == "density"){
            plot1 <-  densityplot(~ tmpvec11 + tmpvec12 ,xlab="Vmax",
                                  plot.points=FALSE,col=cols[1:2],...)
            plot2 <-  densityplot(~ tmpvec21 + tmpvec22 ,xlab="alpha",
                                  plot.points=FALSE,col=cols[1:2],...)
            if(op.level == 2){
                plot3 <-  densityplot(~ tmpvec31 + tmpvec32 ,xlab="rd",
                                      plot.points=FALSE,col=cols[1:2],...)
            }
            if(op.level == 1){
                print(plot1,position=c(0,0,0.5,1),more=TRUE)
                print(plot2,position=c(0.5,0,1,1))
            }else{
                print(plot1,position=c(0,0,0.3333,1),more=TRUE)
                print(plot2,position=c(0.3333,0,0.6666,1),more=TRUE)
                print(plot3,position=c(0.6666,0,1,1))
            }            
        }
    }else
{
    n1 <- x$niter
    n2 <- x2$niter
    n3 <- x3$niter
    maxchainLength <- max(n1,n2,n3)
    tmpvec11 <- x$resuMC[burnin:n1,1]
    tmpvec12 <- x2$resuMC[burnin:n2,1]
    tmpvec13 <- x3$resuMC[burnin:n3,1]
    tmpvec21 <- x$resuMC[burnin:n1,2]
    tmpvec22 <- x2$resuMC[burnin:n2,2]
    tmpvec23 <- x3$resuMC[burnin:n3,2]
    if(op.level == 2){
        tmpvec31 <- x$resuMC[burnin:n1,3]
        tmpvec32 <- x2$resuMC[burnin:n2,3]
        tmpvec33 <- x3$resuMC[burnin:n3,3]
        ymin3 <- min(c(tmpvec31,tmpvec32,tmpvec33))*0.95
        ymax3 <- max(c(tmpvec31,tmpvec32,tmpvec33))*1.05
    }
    ymin1 <- min(c(tmpvec11,tmpvec12,tmpvec13))*0.95
    ymax1 <- max(c(tmpvec11,tmpvec12,tmpvec13))*1.05
    ymin2 <- min(c(tmpvec21,tmpvec22,tmpvec23))*0.95
    ymax2 <- max(c(tmpvec21,tmpvec22,tmpvec23))*1.05
    if(plot.kind == "trace"){
        plot1 <-  xyplot(tmpvec11 ~ burnin:n1 ,
                         xlim=c(I(burnin-0.05*maxchainLength),I(maxchainLength*1.05)),
                         ylim=c(ymin1,ymax1),
                         xlab = "Iterations", type = "l",
                         ylab = expression(paste("Vmax (",mu,mol," ",m^-2," ",s^-1,")")),
                         panel = function(x,y,...){
                             panel.xyplot(x,y,col=cols[1],...)
                             panel.xyplot(burnin:n2,tmpvec12,col=cols[2],...)
                             panel.xyplot(burnin:n3,tmpvec13,col=cols[3],...)
                         },...)                       
        
        plot2 <-  xyplot(tmpvec21 ~ burnin:n1 ,
                         xlim=c(I(burnin-0.05*maxchainLength),I(maxchainLength*1.05)),
                         ylim=c(ymin2,ymax2),
                         xlab = "Iterations", type = "l",
                         ylab = expression(paste("alpha (",mol," ",m^-1,")")),
                         panel = function(x,y,...){
                             panel.xyplot(x,y,col=cols[1],...)
                             panel.xyplot(burnin:n2,tmpvec22,col=cols[2],...)
                             panel.xyplot(burnin:n3,tmpvec23,col=cols[3],...)
                         },...)
        if(op.level == 2){
            plot3 <-  xyplot(tmpvec31 ~ burnin:n1 ,
                             xlim=c(I(burnin-0.05*maxchainLength),I(maxchainLength*1.05)),
                             ylim=c(ymin3,ymax3),
                             xlab = "Iterations", type = "l",
                             ylab = expression(paste("Rd (",mu,mol," ",m^-2," ",s^-1,")")),
                             panel = function(x,y,...){
                                 panel.xyplot(x,y,col=cols[1],...)
                                 panel.xyplot(burnin:n2,tmpvec32,col=cols[2],...)
                                 panel.xyplot(burnin:n3,tmpvec33,col=cols[3],...)
                             },...)
        }
            if(op.level == 1){
                print(plot1,position=c(0,0,0.5,1),more=TRUE)
                print(plot2,position=c(0.5,0,1,1))
            }else{
                print(plot1,position=c(0,0,0.3333,1),more=TRUE)
                print(plot2,position=c(0.3333,0,0.6666,1),more=TRUE)
                print(plot3,position=c(0.6666,0,1,1))
            }         
    } else
    if(plot.kind == "density"){
        plot1 <-  densityplot(~ tmpvec11 + tmpvec12 + tmpvec13
                              ,xlab="Vmax", plot.points=FALSE,col=cols,...)
        plot2 <-  densityplot(~ tmpvec21 + tmpvec22 + tmpvec23
                              ,xlab="alpha", plot.points=FALSE,col=cols,...)
        if(op.level == 2){
            plot3 <-  densityplot(~ tmpvec31 + tmpvec32 + tmpvec33 ,xlab="rd",
                                  plot.points=FALSE,col=cols,...)
        }
        if(op.level == 1){
            print(plot1,position=c(0,0,0.5,1),more=TRUE)
            print(plot2,position=c(0.5,0,1,1))
        }else{
            print(plot1,position=c(0,0,0.3333,1),more=TRUE)
            print(plot2,position=c(0.3333,0,0.6666,1),more=TRUE)
            print(plot3,position=c(0.6666,0,1,1))
        }            
    }
}
}


predict.MCMCc4photo <- function(x, obs, burnin=1e3, kparm = 0.7, theta = 0.83,
                                beta = 0.93, Catm = 380, b0 = 0.08, b1 = 3,
                                stress = 1, ws = "gs"){

    ## The assumption is that x is an object of class "MCMCc4photo"
    ## obs is assumed to be in the same format as the input for the MCMCc4photo function
    niter <- x$niter

    parm <- x$resuMC[burnin:niter,]
    
    nro <- nrow(obs)
    idx1s <- seq(1, niter * nro, nro)
    resd <- data.frame(assim = rep(NA, (niter - burnin) * nro), qp = NA, niter = NA)

    tmp <- data.frame(assim = NA, qp = rep(NA, nro), iter = NA)

    for(i in 1:(niter - burnin)){

        vmax <- parm[i,1]
        alpha <- parm[i,2]
        rd <- parm[i,3]
        assim <- c4photo(obs[,2], obs[,3], obs[,4],
                         vmax=vmax, alpha=alpha, Rd = rd,
                         kparm = kparm, theta = theta, beta = beta,
                         Catm = Catm, b0 = b0, b1 = b1,
                         stress = stress, ws = ws)$Assim
        tmp$assim <- assim
        tmp$qp <- obs[,2]
        tmp$iter <- i
        idx1 <- idx1s[i]
        idx2 <- idx1 + (nro - 1)
        resd[idx1:idx2,] <- tmp

    }
    resd
}
