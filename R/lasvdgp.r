lasvdgpWorker <- function(X0, design, resp, n0, nn,
                          nfea = min(1000,nrow(design)),
                          nsvd = min(3*nn,nfea), nadd = 1,
                          frac = .9, gstart = 0.001,
                          resvdThres = min(5, nn-n0),
                          every = min(5,nn-n0),
                          maxit=100, verb=0)
{
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
    if(nfea > N || nfea <0) stop("illegal nfea")
    if(n0>nfea || n0 <0) stop("illegal n0")
    if(nn<n0 || nn>nfea) stop("illegal nn")
    if(nsvd<n0 || nsvd > nfea) stop("illegal nsvd")
    if(nadd < 0) stop("illegal nadd")
    if(resvdThres < 0) stop("illegal resvdThres")
    if(every < 0) stop("illegal every")

    out <- .C("lasvdGP_R", as.double(t(X0)), as.double(t(design)),
              as.double(resp), as.integer(M), as.integer(N),
              as.integer(m), as.integer(tlen), as.integer(nn),
              as.integer(n0), as.integer(nfea), as.integer(nsvd),
              as.integer(nadd), as.double(frac), as.double(gstart),
              as.integer(resvdThres), as.integer(every),
              as.integer(maxit), as.integer(verb),
              pmean=double(M*tlen), ps2=double(M*tlen)) #package = lasvdgp

    ret <- list(pmean=matrix(out$pmean,nrow=tlen),ps2=matrix(out$ps2,nrow=tlen))
    return(ret)
}
lasvdgpParallel <- function(X0, design, resp, n0, nn,
                            nfea = min(1000,nrow(design)),
                            nsvd = min(3*nn,nfea), nadd = 1,
                            frac = .9, gstart = 0.001,
                            resvdThres = min(5, nn-n0),
                            every = min(5,nn-n0),
                            maxit=100, verb=0, nthread = 4)
{
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
    if(nfea > N || nfea <0) stop("illegal nfea")
    if(n0>nfea || n0 <0) stop("illegal n0")
    if(nn<n0 || nn>nfea) stop("illegal nn")
    if(nsvd<n0 || nsvd > nfea) stop("illegal nsvd")
    if(nadd < 0) stop("illegal nadd")
    if(resvdThres < 0) stop("illegal resvdThres")
    if(every < 0) stop("illegal every")

    blocks <- blockPar(M,nthread)
    X0par <- lapply(blocks,getRowBlock,X0)
    cl <- parallel::makeCluster(nthread)
    ret <- tryCatch(parallel::parLapply(cl,X0par,lasvdgpWorker,design,
                                        resp,n0,nn,nfea,nsvd,nadd,frac,
                                        gstart,resvdThres,every,maxit,verb),
                    finally=parallel::stopCluster(cl))
    pmean <- matrix(unlist(sapply(ret,`[`,"pmean")),nrow=tlen)
    ps2 <- matrix(unlist(sapply(ret,`[`,"ps2")),nrow=tlen)
    ret <- list(pmean=pmean,ps2=ps2)
    return(ret)
}

lasvdgpCall <- function(X0, design, resp, n0, nn,
                        nfea = min(1000,nrow(design)),
                        nsvd = min(3*nn,nfea), nadd = 1,
                        frac = .9, gstart = 0.001,
                        resvdThres = min(5, nn-n0),
                        every = min(5,nn-n0),numstarts=5,
                        maxit=100, verb=0)
{
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
    if(nfea > N || nfea <0) stop("illegal nfea")
    if(n0>nfea || n0 <0) stop("illegal n0")
    if(nn<n0 || nn>nfea) stop("illegal nn")
    if(nsvd<n0 || nsvd > nfea) stop("illegal nsvd")
    if(nadd < 0) stop("illegal nadd")
    if(resvdThres < 0) stop("illegal resvdThres")
    if(every < 0) stop("illegal every")
    ret <- .Call("lasvdgp_Call",t(X0), t(design), resp,
                 as.integer(nn), as.integer(n0),
                 as.integer(nfea), as.integer(nsvd),
                 as.integer(nadd), as.double(frac),
                 as.double(gstart), as.integer(resvdThres),
                 as.integer(every), as.integer(numstarts),
                 as.integer(maxit),
                 as.integer(verb),  PACKAGE="lasvdgp")
    names(ret) <- c("pmean", "ps2mode","ps2mean","ress2mode","ress2mean","range")
    return(ret)
}
lasvdgpCallParal <- function(X0, design, resp, n0, nn,
                             nfea = min(1000,nrow(design)),
                             nsvd = min(3*nn,nfea), nadd = 1,
                             frac = .9, gstart = 0.001,
                             resvdThres = min(5, nn-n0),
                             every = min(5,nn-n0), numstarts=5,
                             maxit=100, verb=0,
                             nthread = 4)
{
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
    if(nfea > N || nfea <0) stop("illegal nfea")
    if(n0>nfea || n0 <0) stop("illegal n0")
    if(nn<n0 || nn>nfea) stop("illegal nn")
    if(nsvd<n0 || nsvd > nfea) stop("illegal nsvd")
    if(nadd < 0) stop("illegal nadd")
    if(resvdThres < 0) stop("illegal resvdThres")
    if(every < 0) stop("illegal every")
    blocks <- blockPar(M,nthread)
    X0par <- lapply(blocks,getRowBlock,X0)
    cl <- parallel::makeCluster(nthread)
    ret <- tryCatch(parallel::parLapply(cl,X0par,lasvdgpCall,design,
                                        resp,n0,nn,nfea,nsvd,nadd,frac,
                                        gstart,resvdThres,every,numstarts,maxit,verb),
                    finally=parallel::stopCluster(cl))
    pmean <- matrix(unlist(sapply(ret,`[`,"pmean")),nrow=tlen)
    ps2mode <- matrix(unlist(sapply(ret,`[`,"ps2mode")),nrow=tlen)
    ps2mean <- matrix(unlist(sapply(ret,`[`,"ps2mean")),nrow=tlen)
    ress2mode <- unlist(sapply(ret,`[`,"ress2mode"))
    ress2mean <- unlist(sapply(ret,`[`,"ress2mean"))
    range <- do.call(c,sapply(ret,`[`,"range"))
    names(ress2mode) <- NULL
    names(ress2mean) <- NULL
    names(range) <- NULL
    out <- list(pmean=pmean,ps2mode=ps2mode,ps2mean=ps2mean,
                ress2mode=ress2mode,ress2mean=ress2mean,range=range)
    return(out)
}
lasvdgpms <- function(X0, design, resp, n0, nn,
                      nfea = min(1000,nrow(design)),
                      nsvd = min(3*nn,nfea), nadd = 1,
                      frac = .9, gstart = 0.001,
                      resvdThres = min(5, nn-n0),
                      every = min(5,nn-n0),
                      nstarts = 5,
                      maxit=100, verb=0)
{
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
    if(nfea > N || nfea <0) stop("illegal nfea")
    if(n0>nfea || n0 <0) stop("illegal n0")
    if(nn<n0 || nn>nfea) stop("illegal nn")
    if(nsvd<n0 || nsvd > nfea) stop("illegal nsvd")
    if(nadd < 0) stop("illegal nadd")
    if(resvdThres < 0) stop("illegal resvdThres")
    if(every < 0) stop("illegal every")
    if(nstarts < 0) stop("illegal nstarts")
    out <- .C("lasvdGPms_R", as.double(t(X0)), as.double(t(design)),
              as.double(resp), as.integer(M), as.integer(N),
              as.integer(m), as.integer(tlen), as.integer(nn),
              as.integer(n0), as.integer(nfea), as.integer(nsvd),
              as.integer(nadd), as.double(frac), as.double(gstart),
              as.integer(resvdThres), as.integer(every), as.integer(nstarts),
              as.integer(maxit), as.integer(verb),
              pmean=double(M*tlen), ps2=double(M*tlen)) #package = lasvdgp
    ret <- list(pmean=matrix(out$pmean,nrow=tlen),ps2=matrix(out$ps2,nrow=tlen))
    return(ret)
}
lasvdgpmsParal <- function(X0, design, resp, n0, nn,
                           nfea = min(1000,nrow(design)),
                           nsvd = min(3*nn,nfea), nadd = 1,
                           frac = .9, gstart = 0.001,
                           resvdThres = min(5, nn-n0),
                           every = min(5,nn-n0),
                           nstarts = 5,
                           maxit=100, verb=0,
                           nthread = 4)
{
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
    if(nfea > N || nfea <0) stop("illegal nfea")
    if(n0>nfea || n0 <0) stop("illegal n0")
    if(nn<n0 || nn>nfea) stop("illegal nn")
    if(nsvd<n0 || nsvd > nfea) stop("illegal nsvd")
    if(nadd < 0) stop("illegal nadd")
    if(resvdThres < 0) stop("illegal resvdThres")
    if(every < 0) stop("illegal every")
    if(nstarts < 0) stop("illegal nstarts")
    blocks <- blockPar(M,nthread)
    X0par <- lapply(blocks,getRowBlock,X0)
    cl <- parallel::makeCluster(nthread)
    ret <- tryCatch(parallel::parLapply(cl,X0par,lasvdgpms,design,
                                        resp,n0,nn,nfea,nsvd,nadd,frac,
                                        gstart,resvdThres,every,nstarts,maxit,verb),
                    finally=parallel::stopCluster(cl))
    pmean <- matrix(unlist(sapply(ret,`[`,"pmean")),nrow=tlen)
    ps2 <- matrix(unlist(sapply(ret,`[`,"ps2")),nrow=tlen)
    out <- list(pmean=pmean,ps2=ps2)
    return(out)
}
