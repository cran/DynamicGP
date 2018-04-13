lasvdGP <- function(design, resp, X0=design, n0=10, nn=20,
                    nfea = min(1000,nrow(design)),
                    nsvd = nn, nadd = 1, frac = .9, gstart = 0.0001,
                    resvdThres = min(5, nn-n0), every = min(5,nn-n0),
                    nstarts = 5,centralize=FALSE, maxit=100, verb=0,
                    nthread = 4, clutype="PSOCK")
{
    if(.Platform$r_arch=="i386")
    {
        cat("the current version does not support 32bit architecture, exit\n")
        return(NULL)
    }
    if(!is.matrix(design)) stop("design must be a matrix")
    if(!is.matrix(resp)) stop("resp must be a matrix")
    N <- nrow(design)
    m <- ncol(design)
    if(ncol(resp) != N)
        stop("number of design points and responses are not consistent")
    tlen <- nrow(resp)
    if(!is.matrix(X0) && length(X0) != m)
        stop("illegal form of prediction set")
    if(!is.matrix(X0)) X0 <- matrix(X0,ncol=m)
    if(ncol(X0) != m) stop("dimensions of design and prediction set are not consistent")
    M <- nrow(X0)
    if(nthread > 1)
    {
        if(nstarts > 1)
            ret <- lasvdgpmsParal(X0,design,resp,n0,nn,nfea,nsvd,nadd,
                                  frac,gstart,resvdThres,every,nstarts,
                                  centralize,maxit,verb,nthread,clutype)
        else
            ret <- lasvdgpParal(X0,design,resp,n0,nn,nfea,nsvd,nadd,
                                frac,gstart,resvdThres,every,
                                centralize,maxit,verb,nthread,clutype)
    }
    else
    {
        if(nstarts > 1)
            ret <- lasvdgpms(X0,design,resp,n0,nn,nfea,nsvd,nadd,
                             frac,gstart,resvdThres,every,nstarts,
                             centralize,maxit,verb)
        else
            ret <- lasvdgpWorker(X0,design,resp,n0,nn,nfea,nsvd,nadd,
                                 frac,gstart,resvdThres,every,
                                 centralize,maxit,verb)
    }
    return(ret)
}
