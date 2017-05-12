blockPar <- function(num,nblock)
{
    bsize <- floor(num/nblock)
    starts <- seq(1,(nblock-1)*bsize+1,len=nblock)
    ends <- starts+bsize-1
    ends[nblock]=num
    mat <- rbind(starts,ends)
    ret <- as.list(as.data.frame(mat))
}
getRowBlock <- function(bound,mat)
{
    idx <- seq(bound[1],bound[2])
    ret <- mat[idx,]
}
