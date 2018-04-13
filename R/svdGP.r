gpsepms <- function(resp,design,test,nstarts=5,d=NULL,g=0.0001)
{
    din <- ncol(design)
    d <- darg(d,design)
    g <- garg(g,resp)
    parb <- maximinLHS(nstarts,din+1)
    ranb <- parb[,1:din]
    nugb <- parb[,din+1]
    ldmin <- log(d$min)
    ldmax <- log(d$max)
    ranb <- ldmin + ranb*(ldmax-ldmin)
    ranb <- exp(ranb)
    lgmin <- log(g$min)
    lgmax <- log(g$max)
    nugb <- lgmin + nugb*(lgmax-lgmin)
    nugb <- exp(nugb)
    gpis <- llik <- rep(NA,nstarts)
    for(i in 1:nstarts)
    {
        gpis[i] <- newGPsep(design,resp,ranb[i,],nugb[i],TRUE)
        mle <- jmleGPsep(gpis[i],drange=c(d$min,d$max),
                               grange=c(g$min,g$max),
                               dab=d$ab,gab=g$ab)
        llik[i] <- llikGPsep(gpis[i],dab=d$ab,gab=g$ab)
    }
    optidx <- which.max(llik)
    pred <- predGPsep(gpis[optidx],test,TRUE)
    for(i in 1:nstarts) deleteGPsep(gpis[i])
    return(pred)
}
buildBasis <- function(response,percent=.9,numbas=NULL)
{
    svdm <- svd(response)
    if(is.null(numbas))
    {
        cumd <- cumsum(svdm$d)/sum(svdm$d)
        numbas <- min(which(cumd>=percent))
    }
    basis <- svdm$u[,1:numbas]          #basis matrix TT by numbas
    redd <- svdm$d[1:numbas]            #reduced d vector length numbas
    redv <- svdm$v[,1:numbas]           #reduced v matrix nn by mumbas
    basis <- t(t(basis)*redd)            #coefficient matrix nn by numbas
    coeff <- redv
    if(!is.matrix(coeff))
        coeff <- matrix(coeff,ncol=numbas)
    ret <- list(basis=basis,redd=redd,coeff=coeff,
                numbas=numbas)
}
svdGP <- function(design,resp,X0=design,nstarts=5,d=NULL,gstart=0.0001,
                  frac=.9,centralize=FALSE,nthread=4,clutype="PSOCK")
{
    if(.Platform$r_arch=="i386")
    {
        cat("the current version does not support 32bit architecture, exit\n")
        return(NULL)
    }
	lenresp <- length(resp)
    if(centralize)
    {
        rmean <- apply(resp,1,mean)
        resp <- resp-rmean
    }
    lbasis <- buildBasis(resp,frac)
    basis <- lbasis$basis
    numbas <- lbasis$numbas
    coeff <- lbasis$coeff
    resid <- resp-basis%*%t(coeff)
    varres <- drop(crossprod(as.vector(resid)))
    varres <- varres/(lenresp+2)
    cl <- parallel::makeCluster(nthread,type=clutype)
    ret <- tryCatch(parallel::parApply(cl,coeff,2,gpsepms,
                                       design,X0,nstarts,d,gstart),
                    finally=parallel::stopCluster(cl))
    vmean <- matrix(unlist(sapply(ret,`[`,"mean")),nrow=numbas,byrow=TRUE)
    vsigma2 <- matrix(unlist(sapply(ret,`[`,"s2")),nrow=numbas,byrow=TRUE)
    pmean <- basis%*%vmean
    ps2 <- basis^2%*%vsigma2+varres
    if(centralize) pmean <- pmean+rmean
    ret <- list(pmean=pmean,ps2=ps2)
    return(ret)
}
