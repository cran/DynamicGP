lasvdgpWorker <- function(X0, design, resp, n0, nn,
                          nfea = min(1000,nrow(design)),
                          nsvd = nn, nadd = 1,
                          frac = .9, gstart = 0.001,
                          resvdThres = min(5, nn-n0),
                          every = min(5,nn-n0),centralize=FALSE,
                          maxit=100, verb=0)
{
    N <- nrow(design)
    m <- ncol(design)
    tlen <- nrow(resp)
    M <- nrow(X0)
    if(nfea > N || nfea <0) stop("illegal nfea")
    if(n0>nfea || n0 <0) stop("illegal n0")
    if(nn<n0 || nn>nfea) stop("illegal nn")
    if(nsvd<n0 || nsvd > nfea) stop("illegal nsvd")
    if(nadd < 0) stop("illegal nadd")
    if(resvdThres < 0) stop("illegal resvdThres")
    if(every < 0) stop("illegal every")
    if(centralize)
    {
        rmean <- apply(resp,1,mean)
        resp <- (resp-rmean)
    }
    out <- .C("lasvdGP_R", as.double(t(X0)), as.double(t(design)),
              as.double(resp), as.integer(M), as.integer(N),
              as.integer(m), as.integer(tlen), as.integer(nn),
              as.integer(n0), as.integer(nfea), as.integer(nsvd),
              as.integer(nadd), as.double(frac), as.double(gstart),
              as.integer(resvdThres), as.integer(every),
              as.integer(maxit), as.integer(verb),
              pmean=double(M*tlen), ps2=double(M*tlen)) #package = lasvdgp
    pmean <- out$pmean
    ps2 <- out$ps2
    if(centralize) pmean <- pmean+rmean
    ret <- list(pmean=matrix(pmean,nrow=tlen),ps2=matrix(ps2,nrow=tlen))
    return(ret)
}
lasvdgpParal <- function(X0, design, resp, n0, nn,
                         nfea = min(1000,nrow(design)),
                         nsvd = nn, nadd = 1,
                         frac = .9, gstart = 0.001,
                         resvdThres = min(5, nn-n0),
                         every = min(5,nn-n0),centralize=FALSE,
                         maxit=100, verb=0, nthread = 4, clutype="PSOCK")
{
    N <- nrow(design)
    m <- ncol(design)
    tlen <- nrow(resp)
    M <- nrow(X0)
    if(nfea > N || nfea <0) stop("illegal nfea")
    if(n0>nfea || n0 <0) stop("illegal n0")
    if(nn<n0 || nn>nfea) stop("illegal nn")
    if(nsvd<n0 || nsvd > nfea) stop("illegal nsvd")
    if(nadd < 0) stop("illegal nadd")
    if(resvdThres < 0) stop("illegal resvdThres")
    if(every < 0) stop("illegal every")
    ## note that the normalization is performed here, so we do not need to normalize again in lasvdgpWorker
    if(centralize)
    {
        rmean <- apply(resp,1,mean)
        resp <- resp-rmean
    }
    blocks <- blockPar(M,nthread)
    X0par <- lapply(blocks,getRowBlock,X0)
    cl <- parallel::makeCluster(nthread,type=clutype)
    ret <- tryCatch(parallel::parLapply(cl,X0par,lasvdgpWorker,design,
                                        resp,n0,nn,nfea,nsvd,nadd,frac,
                                        gstart,resvdThres,every,FALSE,maxit,verb),
                    finally=parallel::stopCluster(cl))
    pmean <- matrix(unlist(sapply(ret,`[`,"pmean")),nrow=tlen)
    ps2 <- matrix(unlist(sapply(ret,`[`,"ps2")),nrow=tlen)
    if(centralize) pmean <- pmean+rmean
    ret <- list(pmean=pmean,ps2=ps2)
    return(ret)
}

lasvdgpms <- function(X0, design, resp, n0, nn,
                      nfea = min(1000,nrow(design)),
                      nsvd = nn, nadd = 1,
                      frac = .9, gstart = 0.001,
                      resvdThres = min(5, nn-n0),
                      every = min(5,nn-n0),
                      nstarts = 5,centralize=FALSE,
                      maxit=100, verb=0)
{
    N <- nrow(design)
    m <- ncol(design)
    tlen <- nrow(resp)
    M <- nrow(X0)
    if(nfea > N || nfea <0) stop("illegal nfea")
    if(n0>nfea || n0 <0) stop("illegal n0")
    if(nn<n0 || nn>nfea) stop("illegal nn")
    if(nsvd<n0 || nsvd > nfea) stop("illegal nsvd")
    if(nadd < 0) stop("illegal nadd")
    if(resvdThres < 0) stop("illegal resvdThres")
    if(every < 0) stop("illegal every")
    if(nstarts < 0) stop("illegal nstarts")
    if(centralize)
    {
        rmean <- apply(resp,1,mean)
        resp <- resp-rmean
    }
    out <- .C("lasvdGPms_R", as.double(t(X0)), as.double(t(design)),
              as.double(resp), as.integer(M), as.integer(N),
              as.integer(m), as.integer(tlen), as.integer(nn),
              as.integer(n0), as.integer(nfea), as.integer(nsvd),
              as.integer(nadd), as.double(frac), as.double(gstart),
              as.integer(resvdThres), as.integer(every), as.integer(nstarts),
              as.integer(maxit), as.integer(verb),
              pmean=double(M*tlen), ps2=double(M*tlen)) #package = lasvdgp
    pmean <- out$pmean
    ps2 <- out$ps2
    if(centralize) pmean <- pmean+rmean
    ret <- list(pmean=matrix(pmean,nrow=tlen),ps2=matrix(ps2,nrow=tlen))
    return(ret)
}
lasvdgpmsParal <- function(X0, design, resp, n0, nn,
                           nfea = min(1000,nrow(design)),
                           nsvd = nn, nadd = 1,
                           frac = .9, gstart = 0.001,
                           resvdThres = min(5, nn-n0),
                           every = min(5,nn-n0),
                           nstarts = 5,centralize=FALSE,
                           maxit=100, verb=0,
                           nthread = 4, clutype="PSOCK")
{
    N <- nrow(design)
    m <- ncol(design)
    tlen <- nrow(resp)
    M <- nrow(X0)
    if(nfea > N || nfea <0) stop("illegal nfea")
    if(n0>nfea || n0 <0) stop("illegal n0")
    if(nn<n0 || nn>nfea) stop("illegal nn")
    if(nsvd<n0 || nsvd > nfea) stop("illegal nsvd")
    if(nadd < 0) stop("illegal nadd")
    if(resvdThres < 0) stop("illegal resvdThres")
    if(every < 0) stop("illegal every")
    if(nstarts < 0) stop("illegal nstarts")
    ## note that the normalization is performed here, so we do not need to normalize again in lasvdgpms
    if(centralize)
    {
        rmean <- apply(resp,1,mean)
        resp <- resp-rmean
    }
    blocks <- blockPar(M,nthread)
    X0par <- lapply(blocks,getRowBlock,X0)
    cl <- parallel::makeCluster(nthread,type=clutype)
    ret <- tryCatch(parallel::parLapply(cl,X0par,lasvdgpms,design,
                                        resp,n0,nn,nfea,nsvd,nadd,frac,
                                        gstart,resvdThres,every,nstarts,
                                        FALSE,maxit,verb),
                    finally=parallel::stopCluster(cl))
    pmean <- matrix(unlist(sapply(ret,`[`,"pmean")),nrow=tlen)
    ps2 <- matrix(unlist(sapply(ret,`[`,"ps2")),nrow=tlen)
    if(centralize) pmean <- pmean+rmean
    out <- list(pmean=pmean,ps2=ps2)
    return(out)
}
