### This file contains the function for the dots plot.
### Plot the difference for all sequences based on the first sequence.

plotdots <- function(X, X.class = NULL, Mu = NULL, code.type = .code.type[1],
    diff.only = TRUE, fill = FALSE, label = TRUE, xlim = NULL,
    main = "Dots Plot", xlab = "Sites", ylab = "Sequences", ...){
  if(sum(code.type %in% .code.type) != 1){
    stop("code.type is not found.")
  }

  if(! is.null(X.class)){
    if(length(X.class) != nrow(X)){
      stop("length(X.class) != N.")
    }

    id.X <- order(X.class)
    X <- X[id.X,]
    K <- max(X.class)
    X.nc.cum <- rep(0, K)
    for(k in 1:K) X.nc.cum[k] <- sum(X.class %in% 1:k)
  }

  if(is.null(Mu)){
    Mu <- find.consensus(X, code.type = code.type)
  } else{
    if(length(Mu) != ncol(X)){
      stop("length(Mu) != L.")
    }
  }

  if(diff.only){
    tl.diff <- apply(X, 2, function(x) length(unique(x)))
    X <- X[, tl.diff > 1]
    Mu <- Mu[tl.diff > 1]
  }

  N <- nrow(X)
  L <- ncol(X)

  if(is.null(xlim)) xlim <- c(1, L + 1)
  if(label){
    ylim <- c(N + 3, -1)
  } else{
    ylim <- c(N + 1, -1)
  }

  plot(NULL, NULL, type = "n", xlim = xlim, ylim = ylim,
       main = main, xlab = xlab, ylab = ylab)

  my.col <- c("green3", "blue2", "#CC00CC", "red2", "black")
  flag.diff <- rowSums(t(X) != Mu) > 0

  for(j in 1:L){
    x.left <- j
    y.top <- -1
    rect(x.left, y.top + 2, x.left + 1, y.top,
         col = my.col[Mu[j] + 1], border = NA)
  }
  abline(h = y.top + 2, lty = 3, lwd = 0.5)

  for(i in 1:N){
    for(j in 1:L){
      if(diff.only){
         if(! flag.diff[j]) next
      }
      x.left <- j
      y.top <- i
      if(X[i, j] != Mu[j]){
        rect(x.left, y.top + 1, x.left + 1, y.top,
             col = my.col[X[i, j] + 1], border = NA)
      } else{
        if(fill & flag.diff[j]){
          rect(x.left, y.top + 1, x.left + 1, y.top,
               col = my.col[X[i, j] + 1], border = NA)
        }
      }
    }
  }

  if(label){
    i <- N + 1
    for(j in 1:L){
      x.left <- j
      y.top <- i
      if(flag.diff[j]){
        rect(x.left, y.top + 2, x.left + 1, y.top, col = "orange", border = NA)
      }
    }
    abline(h = y.top, lty = 3, lwd = 0.5)
  }

  if(! is.null(X.class)) abline(h = X.nc.cum + 1, lty = 3, lwd = 0.5)
} # End of plotdots().
