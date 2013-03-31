
grape <- function(nnode, nbatch, state, theta, nspac = 1, blen = 1)
{
    if (length(nnode) != 1)
        stop("nnode not scalar")
    if (length(nbatch) != 1)
        stop("nbatch not scalar")
    if (length(blen) != 1)
        stop("blen not scalar")
    if (length(nspac) != 1)
        stop("nspac not scalar")

    inode <- as.integer(nnode)
    ibatch <- as.integer(nbatch)
    iblen <- as.integer(blen)
    ispac <- as.integer(nspac)
    istate <- state
    storage.mode(istate) <- "integer"
    dtheta <- as.double(theta)

    if (inode != nnode)
        stop("nnode not integer")
    if (ibatch != nbatch)
        stop("nbatch not integer")
    if (iblen != blen)
        stop("blen not integer")
    if (ispac != nspac)
        stop("nspac not integer")
    if (! all(istate == state))
        stop("state not integer")

    if (inode <= 0)
        stop("nnode not positive")
    if (ibatch <= 0)
        stop("nbatch not positive")
    if (iblen <= 0)
        stop("blen not positive")
    if (ispac <= 0)
        stop("nspac not positive")

    if (length(dtheta) != inode)
        stop("theta wrong dimension (not length nnode)")
    if (! is.matrix(istate))
        stop("state not matrix")
    if (! all(dim(istate) == inode))
        stop("state not nnode by nnode")

    if (! all(istate == 0 | istate == 1))
        stop("state not 0-or-1-valued")
    if (! all(diag(istate) == 0))
        stop("state not zero diagonal")
    if (! all(istate == t(istate)))
        stop("state not symmetric")

    dbatch <- matrix(as.double(0), nnode, nbatch)

    out.time <- system.time(
    out <- .C("grape", nnode = inode, nspac = ispac, nbatch = ibatch,
        blen = iblen, state = istate, theta = dtheta, batch = dbatch,
        PACKAGE = "grape")
    )

    out$state[upper.tri(out$state, diag = TRUE)] <- 0
    out$state <- out$state + t(out$state)
    out$batch <- t(out$batch)
    out$time <- out.time
    return(out)
}

