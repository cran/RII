".required" <-
"splines"
"LL" <-
function (para, count, classcount, agecount)
## Log likelihood for linear rate 
{
    MU <- (MU.co.sl %*% para[1:2] %*% t(c(1, exp(para[-(1:2)])))) * 
        n.mat
    if (any(MU < 0)) 
        return(Inf)
    -sum(count * log(MU)) + sum(MU)
}
"LL.grad" <-
function (para, count, classcount, agecount) 
## Gradient of log likelihood for linear rate
{
    c(crossprod(MU.co.sl, n.mat %*% c(1, exp(para[-(1:2)])) - 
        classcount/(MU.co.sl %*% para[1:2])), crossprod(n.mat[, 
        -1], MU.co.sl %*% para[1:2]) * exp(para[-(1:2)]) - agecount)
}
"LL.more" <-
function (para, count, classcount, agecount)
## Log likelihood for linear rate; also returns fitted counts 
{
    MU <- (MU.co.sl %*% para[1:2] %*% t(c(1, exp(para[-(1:2)])))) * 
        n.mat
    z <- -sum(count * log(MU)) + sum(MU)
    list(val = z, MU = MU)
}
"myns" <-
function (x, df = NULL, knots = NULL, intercept = F, Boundary.knots = range(x), 
    derivs = 0) 
{
    nx <- names(x)
    x <- as.vector(x)
    nax <- is.na(x)
    if (nas <- any(nax)) 
        x <- x[!nax]
    if (!missing(Boundary.knots)) {
        if (length(Boundary.knots) == 1 && Boundary.knots == 
            T) {
            Boundary.knots <- range(knots)
            knots <- sort(knots)[-c(1, length(knots))]
        }
        Boundary.knots <- range(Boundary.knots)
        outside <- (ol <- x < Boundary.knots[1]) | (or <- x > 
            Boundary.knots[2])
    }
    else outside <- rep(F, length = length(x))
    if (!missing(df) && missing(knots)) {
        nIknots <- df - 1 - intercept
        if (nIknots < 0) {
            nIknots <- 0
            warning(paste("df was too small; have used ", 1 + 
                intercept))
        }
        if (nIknots > 0) {
            knots <- seq(from = 0, to = 1, length = nIknots + 
                2)[-c(1, nIknots + 2)]
            knots <- quantile(x[!outside], knots)
        }
        else knots <- NULL
    }
    Aknots <- sort(c(rep(Boundary.knots, 4), knots))
    if (any(outside)) {
        basis <- array(0, c(length(x), length(knots) + 4))
        if (any(ol)) {
            k.pivot <- Boundary.knots[1]
            xl <- cbind(1, x[ol] - k.pivot)
            tt <- spline.des(Aknots, rep(k.pivot, 2), 4, c(derivs, 
                2))$design
            basis[ol, ] <- xl %*% tt
        }
        if (any(or)) {
            k.pivot <- Boundary.knots[2]
            xr <- cbind(1, x[or] - k.pivot)
            tt <- spline.des(Aknots, rep(k.pivot, 2), 4, c(derivs, 
                2))$design
            basis[or, ] <- xr %*% tt
        }
        if (any(inside <- !outside)) 
            basis[inside, ] <- spline.des(Aknots, x[inside], 
                4, rep(derivs, length(x[inside])))$design
    }
    else basis <- spline.des(Aknots, x, 4, rep(derivs, length(x)))$design
    const <- spline.des(Aknots, Boundary.knots, 4, c(2, 2))$design
    if (!intercept) {
        const <- const[, -1, drop = F]
        basis <- basis[, -1, drop = F]
    }
    qr.const <- qr(t(const))
    basis <- as.matrix((t(qr.qty(qr.const, t(basis))))[, -(1:2)])
    n.col <- ncol(basis)
    if (nas) {
        nmat <- matrix(NA, length(nax), n.col)
        nmat[!nax, ] <- basis
        basis <- nmat
    }
    dimnames(basis) <- list(nx, 1:n.col)
    a <- list(degree = 3, knots = knots, Boundary.knots = Boundary.knots, 
        intercept = intercept, class = c("ns", "basis"))
    attributes(basis) <- c(attributes(basis), a)
    basis
}
"ncs.basis" <-
function (knots) 
{
    nknot <- length(knots)
    basis.mat <- matrix(ns(knots, knots = knots[2:(nknot - 1)], 
        intercept = TRUE, Boundary.knots = c(0, 1)), nrow = nknot)
    derivs.mat <- matrix(myns(knots, knots = knots[2:(nknot - 
        1)], intercept = TRUE, Boundary.knots = c(0, 1), derivs = 2), 
        nrow = nknot)
    A0 <- A1 <- B0 <- B1 <- B2 <- B3 <- matrix(0, nknot - 1, 
        nknot)
    for (i in 1:(nknot - 1)) {
        h <- knots[i + 1] - knots[i]
        for (j in 1:nknot) {
            tL <- knots[i]
            tR <- knots[i + 1]
            fL <- basis.mat[i, j]
            fR <- basis.mat[i + 1, j]
            dL <- derivs.mat[i, j]
            dR <- derivs.mat[i + 1, j]
            A0[i, j] <- (tR * dL - tL * dR)/h
            A1[i, j] <- (dR - dL)/h
            B0[i, j] <- (dL * tR^3 - dR * tL^3)/(6 * h) + (tR * 
                fL - tL * fR)/h + h * (dR * tL - dL * tR)/6
            B1[i, j] <- (dR * tL^2 - dL * tR^2)/(2 * h) + (fR - 
                fL)/h + h * (dL - dR)/6
            B2[i, j] <- (dL * tR - dR * tL)/(2 * h)
            B3[i, j] <- (dR - dL)/(6 * h)
        }
    }
    list(A0 = A0, A1 = A1, B0 = B0, B1 = B1, B2 = B2, B3 = B3)
}
"PLL" <-
function (para, count, classcount, agecount, L)
## Penalized log likelihood 
{
    MU <- (MU.co %*% para[1:p] %*% t(c(1, exp(para[-(1:p)])))) * 
        n.mat
    if (any(MU < 0)) 
        return(Inf)
    penalty <- pen.co * (para[1:p] %*% t(para[1:p]))
    -sum(count * log(MU)) + sum(MU) + exp(L) * sum(penalty)
}
"PLL.ex" <-
function (para, count, classcount, L, ex)
## Penalized log likelihood excluding class `ex' 
{
    MU <- (MU.co %*% para[1:p] %*% t(c(1, exp(para[-(1:p)])))) * 
        n.mat
    if (any(MU < 0)) 
        return(Inf)
    penalty <- pen.co * (para[1:p] %*% t(para[1:p]))
    -sum(count[-ex, ] * log(MU[-ex, ])) + sum(MU[-ex, ]) + exp(L) * 
        sum(penalty)
}
"PLL.ex.grad" <-
function (para, count, classcount, L, ex)
## Gradient of penalized log likelihood excluding class `ex' 
{
    c(crossprod(MU.co[-ex, ], n.mat[-ex, , drop = FALSE] %*% 
        c(1, exp(para[-(1:p)])) - classcount[-ex]/(MU.co[-ex, 
        ] %*% para[1:p])), crossprod(n.mat[-ex, -1], MU.co[-ex, 
        ] %*% para[1:p]) * exp(para[-(1:p)]) - colSums(count[-ex, 
        , drop = FALSE])[-1]) + c(exp(L) * ((pen.co %*% para[1:p]) + 
        crossprod(pen.co, para[1:p])), rep(0, nage - 1))
}
"PLL.ex.more" <-
function (para, count, classcount, L, ex)
## Penalized log likelihood excluding class `ex'; also returns fitted counts 
{
    MU <- (MU.co %*% para[1:p] %*% t(c(1, exp(para[-(1:p)])))) * 
        n.mat
    penalty <- pen.co * (para[1:p] %*% t(para[1:p]))
    z <- -sum(count[-ex, ] * log(MU[-ex, ])) + sum(MU[-ex, ]) + 
        exp(L) * sum(penalty)
    list(val = z, MU = MU)
}
"PLL.grad" <-
function (para, count, classcount, agecount, L)
## Gradient of penalized log likelihood 
{
    c(crossprod(MU.co, n.mat %*% c(1, exp(para[-(1:p)])) - classcount/(MU.co %*% 
        para[1:p])), crossprod(n.mat[, -1], MU.co %*% para[1:p]) * 
        exp(para[-(1:p)]) - agecount) + c(exp(L) * ((pen.co %*% 
        para[1:p]) + crossprod(pen.co, para[1:p])), rep(0, nage - 
        1))
}
"PLL.more" <-
function (para, count, classcount, agecount, L)
## Penalized log likelihood; also returns fitted counts 
{
    MU <- (MU.co %*% para[1:p] %*% t(c(1, exp(para[-(1:p)])))) * 
        n.mat
    penalty <- pen.co * (para[1:p] %*% t(para[1:p]))
    z <- -sum(count * log(MU)) + sum(MU) + exp(L) * sum(penalty)
    list(val = z, MU = MU, penalty = penalty)
}
"plot.RII" <-
function (x, group = 1, ylim = NULL, ...) 
{
    z <- x
    nknot <- dim(z$pop)[1]
    knots <- seq(0, 1, len = nknot)
    bounds <- c(0, cumsum(rowSums(z$pop))/sum(z$pop))
    x <- sort(c(bounds[1], rep(bounds[2:nknot], 2), bounds[nknot + 
        1]))
    rate <- as.vector(rbind(z$count[, group]/z$pop[, group], 
        z$count[, group]/z$pop[, group]))
    basis <- matrix(ns(seq(0, 1, len = 101), knots = knots[2:(nknot - 
        1)], intercept = TRUE, Boundary.knots = c(0, 1)), nrow = 101)
    points.spline <- (basis %*% z$par) * exp(z$group.effects[group])
    if (is.null(ylim)) 
        ylim <- range(c(rate, points.spline))
    plot(x, rate, type = "l", ylim = ylim, ...)
    lines(seq(0, 1, len = 101), points.spline)
}
"RII" <-
function (count, pop, loglambda = NULL, grid = NULL, se = FALSE, 
    B = NULL, alpha = 0.025, returnboot = FALSE) 
{
    name <- dimnames(count)
    count.orig <- count
    pop.orig <- pop
    count <- matrix(count, dim(pop)[1])
    pop <- matrix(pop, dim(pop)[1])
    setup <- RII.setup(pop, dim(pop)[1])
    attach(setup)
    if (se == TRUE) {
        if (B < 1) 
            stop("se is TRUE requires B >= 1")
        boot.data <- array(rpois(nclass * nage * B, count), dim = c(nclass, 
            nage, B))
        boot.rep <- vector(len = B)
        if (is.null(loglambda)) {
            if (is.null(grid)) 
                stop("loglambda = NULL requires grid != NULL")
            loglambda <- RII.loglambda(count, pop, grid, TRUE, boot.data, B)
            ans <- RII.optim(count, pop, loglambda[1], TRUE)
            for (i in 1:B) boot.rep[i] <- RII.optim(as.matrix(boot.data[, 
                , i]), pop, loglambda[i + 1], FALSE)
        }
        else {
            ans <- RII.optim(count, pop, loglambda, TRUE)
            for (i in 1:B) boot.rep[i] <- RII.optim(as.matrix(boot.data[, 
                , i]), pop, loglambda, FALSE)
        }
        ans$se <- sd(log(boot.rep))
        ans$alpha <- alpha
        if(B*alpha >= 1){
            if(B*alpha - floor(B*alpha) < 1e-16)
                ans$ci <- sort(boot.rep)[c(B*alpha,B*(1-alpha))]
            else {
                k <- floor((B + 1) * alpha)
                ans$ci <- sort(boot.rep)[c(k, B - 1 - k)]
            }
            names(ans$ci) <- c(paste(100*alpha,"%"), paste(100*(1-alpha),"%"))
        }
        else cat("Need B*alpha >= 1 to produce a confidence interval")
        if(returnboot == TRUE){
            ans$boot.data <- boot.data
            ans$boot.rep <- boot.rep
            ans$boot.lambda <- loglambda
        }
    }
    else {
        if (is.null(loglambda)) {
            if (is.null(grid)) 
                stop("loglambda = NULL requires grid != NULL")
            loglambda <- RII.loglambda(count, pop, grid = grid)
        }
        ans <- RII.optim(count, pop, loglambda, TRUE)
    }
    detach(setup)
    names(ans$group.effects) <- colnames(count.orig)
    dimnames(ans$count) <- dimnames(ans$pop) <- name
    dimnames(ans$expected) <- dimnames(ans$residuals) <- name
    class(ans) <- "RII"
    ans
}
"RII.CVplot" <-
function (count, pop, loglambda, ...) 
{
    classcount <- rowSums(count)
    setup <- RII.setup(pop, dim(pop)[1] - 1)
    attach(setup)
    temp <- lm(classcount/classpop ~ mids, weights = classpop)$coefficients
    start.points <- temp[1] + temp[2] * mids
    para.start <- c(as.vector(lm(start.points ~ -1 + basis.start)$coefficients), 
        rep(0, nage - 1))
    score <- sapply(loglambda, scorefunc, count = count, classcount = classcount, 
        para.start = para.start)
    detach(setup)
    plot(loglambda, score, ...)
}
"RII.loglambda" <-
function (count, pop, grid, resamp = FALSE, resamp.count, 
    B)
## Returns value of lamdba that minimizes CV score 
{
    setup <- RII.setup(pop, dim(pop)[1] - 1)
    attach(setup)
    RII.CV <- function(count, pop, grid) {
        classcount <- rowSums(count)
        temp <- lm(classcount/classpop ~ mids, weights = classpop)$coefficients
        start.points <- temp[1] + temp[2] * mids
        para.start <- c(as.vector(lm(start.points ~ -1 + basis.start)$coefficients), 
            rep(0, nage - 1))
        score <- sapply(grid, scorefunc, count = count, classcount = classcount, 
            para.start = para.start)
        start <- grid[which.min(score)]
        if (start == min(grid)) 
            return(-Inf)
        if (start == max(grid)) 
            return(Inf)
        optim(start, scorefunc, method = "BFGS", control = list(maxit = 500, 
            fnscale = 1, reltol = 1e-16), count = count, classcount = classcount, 
            para.start = para.start)$par
    }
    loglambda <- RII.CV(count, pop, grid)
    if (resamp == TRUE) {
        for (i in 1:B) loglambda[i + 1] <- RII.CV(resamp.count[, 
            , i], pop, grid = grid)
    }
    detach(setup)
    loglambda
}
"RII.optim" <-
function (count, pop, loglambda, delta = FALSE)
## Returns the value of the RII 
{
    classcount <- rowSums(count)
    agecount <- colSums(count)[-1]
    if (loglambda == Inf) {
        para.start <- c(as.vector(lm(classcount/classpop ~ mids, 
            weights = classpop)$coefficients), rep(0, nage - 
            1))
        para.optim <- optim(para.start, LL, gr = LL.grad, method = "BFGS", 
            control = list(maxit = 500, fnscale = 1, reltol = 1e-16), 
            count = count, classcount = classcount, agecount = agecount)$par
        RII <- para.optim[1]/(para.optim[1] + para.optim[2])
        if (delta == TRUE) {
            temp <- LL.more(para.optim, count = count, classcount = classcount, 
                agecount = agecount)
            maxval <- -temp$val
            expected <- temp$MU
            dev.resids <- sign(count - expected) * sqrt(2 * (count * 
                log(count/expected) - (count - expected)))
            npara <- 1 + nage
            info.mat <- matrix(0, nrow = npara, ncol = npara)
            for (r in 1:2) {
                for (s in 1:2) {
                  info.mat[r, s] <- sum(classcount * MU.co.sl[, 
                    r] * MU.co.sl[, s]/(MU.co.sl %*% para.optim[1:2])^2)
                }
            }
            if (nage > 1) {
                for (r in 3:npara) {
                  info.mat[r, r] <- sum(n.mat[, r - 1] * exp(para.optim[r]) * 
                    (MU.co.sl %*% para.optim[1:2]))
                }
                for (r in 1:2) {
                  for (s in 3:npara) {
                    info.mat[r, s] <- info.mat[s, r] <- sum(n.mat[, 
                      s - 1] * exp(para.optim[s]) * MU.co.sl[, 
                      r])
                  }
                }
            }
            inv.mat <- chol2inv(chol(info.mat))[1:2, 1:2]
            g.dash <- c(para.optim[2], -para.optim[1])/(para.optim[1] + 
                para.optim[2])^2
            variance <- sum((g.dash %*% t(g.dash)) * inv.mat)
            g.log.dash <- c(1/para.optim[1] - 1/(para.optim[1] + 
                para.optim[2]), -1/(para.optim[1] + para.optim[2]))
            variance.log <- sum((g.log.dash %*% t(g.log.dash)) * 
                inv.mat)
        }
    }
    else {
        temp <- lm(classcount/classpop ~ mids, weights = classpop)$coefficients
        start.points <- temp[1] + temp[2] * mids
        para.start <- c(as.vector(lm(start.points ~ -1 + basis.start)$coefficients), 
            rep(0, nage - 1))
        para.optim <- optim(para.start, PLL, gr = PLL.grad, method = "BFGS", 
            control = list(maxit = 500, fnscale = 1, reltol = 1e-16), 
            count = count, classcount = classcount, agecount = agecount, 
            L = loglambda)$par
        points <- basis.end %*% para.optim[1:p]
        RII <- points[1]/points[2]
        if (delta == TRUE) {
            temp <- PLL.more(para.optim, count = count, classcount = classcount, 
                agecount = agecount, L = loglambda)
            maxval <- -temp$val
            expected <- temp$MU
            dev.resids <- sign(count - expected) * sqrt(2 * (count * 
                log(count/expected) - (count - expected)))
            npara <- p + nage - 1
            info.mat <- matrix(0, nrow = npara, ncol = npara)
            for (r in 1:p) {
                for (s in 1:p) {
                  info.mat[r, s] <- sum(classcount * MU.co[, 
                    r] * MU.co[, s]/(MU.co %*% para.optim[1:p])^2)
                }
            }
            if (nage > 1) {
                for (r in (p + 1):npara) {
                  info.mat[r, r] <- sum(n.mat[, r - p + 1] * 
                    exp(para.optim[r]) * (MU.co %*% para.optim[1:p]))
                }
                for (r in 1:p) {
                  for (s in (p + 1):npara) {
                    info.mat[r, s] <- info.mat[s, r] <- sum(n.mat[, 
                      s - p + 1] * exp(para.optim[s]) * MU.co[, 
                      r])
                  }
                }
            }
            inv.mat <- chol2inv(chol(info.mat))[1:p, 1:p]
            g.dash <- (points[2] * basis.end[1, ] - points[1] * 
                basis.end[2, ])/(points[2]^2)
            variance <- sum((g.dash %*% t(g.dash)) * inv.mat)
            g.log.dash <- basis.end[1, ]/points[1] - basis.end[2, 
                ]/points[2]
            variance.log <- sum((g.log.dash %*% t(g.log.dash)) * 
                inv.mat)
        }
    }
    if (delta == FALSE) 
        return(RII)
    npar <- npara - nage + 1
    par = para.optim[1:npar]
    if(npar == 2)
        names(par) <- c("intercept","gradient")
    group.effects = c(0, para.optim[-(1:npar)])
    model <- list(count = count, pop = pop, loglambda = loglambda, 
        par = par, group.effects = group.effects, 
        maxval = maxval, expected = expected, residuals = dev.resids, 
        var = variance, var.log = variance.log, RII = RII)
    model
}
"RII.setup" <-
function (pop, nknot) 
## Does the spadework for RII functions 
{
    knots <- seq(0, 1, len = nknot)
    classpop <- rowSums(pop)
    nclass <- dim(pop)[1]
    nage <- dim(pop)[2]
    expected <- matrix(nrow = nclass, ncol = nage)
    p <- nknot
    nintervals <- nknot - 1
    bounds <- c(0, cumsum(classpop)/sum(classpop))
    mids <- (bounds[-1] + bounds[1:nclass])/2
    breaks <- sort(union(bounds, knots))
    nbreak <- length(breaks)
    breakpos <- vector(len = nbreak - 1)
    for (i in 2:nbreak) breakpos[i - 1] <- max(which(knots < 
        breaks[i]))
    intpos <- vector(len = nbreak - 1)
    for (i in 2:nbreak) intpos[i - 1] <- max(which(bounds < breaks[i]))
    basis.start <- matrix(ns(mids, knots = knots[2:(nknot - 1)], 
        intercept = TRUE, Boundary.knots = c(0, 1)), nrow = nclass)
    basis.end <- matrix(ns(c(0, 1), knots = knots[2:(nknot - 
        1)], intercept = TRUE, Boundary.knots = c(0, 1)), nrow = 2)
    basis <- ncs.basis(knots)
    attach(basis)
    info.ints <- B0[breakpos, , drop = FALSE] * diff(breaks) + 
        B1[breakpos, , drop = FALSE] * diff(breaks^2)/2 + B2[breakpos, 
        , drop = FALSE] * diff(breaks^3)/3 + B3[breakpos, , drop = FALSE] * 
        diff(breaks^4)/4
    info.classints <- matrix(nrow = nclass, ncol = p)
    for (i in 1:nclass) info.classints[i, ] <- colSums(info.ints[which(intpos == 
        i), , drop = FALSE])
    MU.co <- info.classints * sum(pop)
    MU.co.sl <- sum(pop) * cbind(diff(bounds), diff(bounds^2)/2)
    pen.co <- crossprod(A0, A0 * diff(knots)) + crossprod(A0, 
        A1 * diff(knots^2)) + crossprod(A1, A1 * diff(knots^3)/3)
    n.mat <- pop/classpop
    detach(basis)
    list(classpop = classpop, nclass = nclass, nage = nage, nknot = nknot, 
        p = p, expected = expected, mids = mids, basis.start = basis.start, 
        basis.end = basis.end, MU.co = MU.co, MU.co.sl = MU.co.sl, 
        n.mat = n.mat, pen.co = pen.co)
}
"scorefunc" <-
function (loglambda, count, classcount, para.start)
## Returns the CV score obtained using `loglambda' 
{
    for (i in 1:nclass) {
        para.optim <- optim(para.start, PLL.ex, gr = PLL.ex.grad, 
            method = "BFGS", control = list(maxit = 500, fnscale = 1, 
                reltol = 1e-16), count = count, classcount = classcount, 
            L = loglambda, ex = i)$par
        expected[i,] <- PLL.ex.more(para.optim, count = count, 
            L = loglambda, ex = i)$MU[i, ]
    }
    2 * sum(count * log(count/expected) - (count - 
        expected))
}
print.RII <-
function (x, digits = max(3, getOption("digits") - 3),
          na.print = "", 
          ...)
{
    cat("RII estimate:", format(round(x$RII, digits), nsmall = 2), "\n")
    if (length(x$group.effects) > 1){
        cat("\nGroup effects:\n")
        print(round(x$group.effects, digits))
        }
    invisible(x)
}
summary.RII <-
function (object, digits = max(3, .Options$digits - 3), ...)
{
    if(!is.null(object$se)){
        summ <- c(object$RII, object$ci, log(object$RII), object$se)
        names(summ) <- c("RII", names(object$ci), "log(RII)", "se")
        object$summary <- summ
        object$digits <- digits
        class(object) <- "summary.RII"
    }
    object
}
print.summary.RII <-
function (x, digits = x$digits, ...)
{
    cat("RII estimate:\n")
    print(round(x$summary,digits))
    if (length(x$group.effects) > 1){
        cat("\nGroup effects:\n")
        print(round(x$group.effects, digits))
        }
    invisible(x)
}
