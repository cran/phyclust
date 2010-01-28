## modified on Sep 22, 2006. changes in "getcut.fun"
## add "base=2 for calculating "info"

############################################
## FUNCTIONS needed for Hap-RD
############################################

getcut.fun<-function(pp.org,nn=2*nhap, plot=0){
      pp     <-rev(sort(pp.org))
      ct     <-round(pp*nn)
      dimen  <-log((1:length(pp)), base=2 )/ct
      info   <-cumsum(pp*log(1/pp,base=2))
      netinfo<-info-dimen
      cutpos <-netinfo==max(netinfo)
      if(plot==1){
        plot( netinfo, pp,type="b",cex=0.5)
        barplot(pp.org)
        abline(h=pp[cutpos])
      }
      return( (1:length(pp))[cutpos])
    }

# source("code.getPIstar.getBigMatB")	### Remarked by Wei-Chen Chen 2009-09-16

    final.BigMatB.matC.fun<-function(mmle.cscn, hapreserv, loci, digit){
      ## HERE mmle.cs, mmle.cn MUST NOT CONTAIN "0" CATEGORIES!!!
      ## hapreserv <-sort(unique(c(names(mmle.cs[mmle.cs>=cut.cs]),names(mmle.cn[mmle.cn>=cut.cn]))))
      ini.cscn     <-get.ini.possHap.hapReserv.fun(mmle.cscn, loci, digit, hapreserv)
      mut.mat.cscn <-ini.cscn$prepsi
      subPI.cscn   <-get.subPI.fun(mmle.cscn,  mut.mat=mut.mat.cscn)
      ## get BB matrix
      jnk.cscn     <-get.preBB.skip.fun(subPI.cscn, loci, digit, mmle.cscn)
      preBB.cscn   <-jnk.cscn$preBB.lst
      indexBB.cscn <-jnk.cscn$indexBB
      medBB.cscn   <-get.medBB.skip.fun(subPI.lst=subPI.cscn, preBB.lst=preBB.cscn, indexBB.lst=indexBB.cscn,loci, digit, mmle.cscn)
      fullBB.cscn  <-get.fullBB.fun(medBB.lst=medBB.cscn,  subPI.lst=subPI.cscn)
      BB.cscn      <-get.BB.fun(fullBB.cscn)
      CC.cscn      <-lapply(BB.cscn, function(mat){1*(mat>0)})
      BigMatB.cscn <-getBigMatB.fun(BB.cscn, subPI.cscn)
      return(list(## fullBB.cscn=fullBB.cscn,
                  BB.cscn=BB.cscn,
                  CC.cscn =CC.cscn,
                  BigMatB.cscn =BigMatB.cscn ,
                  subPI.cscn =  subPI.cscn
                  ))
  }


# source("functions.r")	### Remarked by Wei-Chen Chen 2009-09-16
## call  all fun related to calcuating deri.BigB and I13.I33i.I13t
## description is in /home/jytzeng/Research/Hap-RD/Score/MatB.cluster/Fang/generalform.ps



get.matD.matEI.fun<-function(yy=y, mmu=mu, aa=a,
                            XRD.post=xRD.post, XFD.post=x.post, XRD=xRD, XFD = x,
                            post.prob=post, V=v, Nreps=nreps,
                            CC.lst=CC, subPI.lst=subPI, BB.lst=BB,
                            pFD=p.cscn.FD, xx.adj=x.adj
                            ){
  ## this use a hybird I matrix: that is, to use obs I_G for all nonzefo element of E(I_G). So Iaz.2=Iaz.3=0,
  ## and Izz.12= Izz.13=Izz.23 =0 
  ## pFD must match with nrow(matB)
  res    <- ((yy - mmu)/aa)^2
  uRD.mtx<-  (yy - mmu)/aa * XRD.post
  RR     <- ncol(XFD)
  ##-----
  ## Daa
  ##-----
  Daa    <- t(uRD.mtx) %*% uRD.mtx  ## sum_i res * E(Xi*|G) * E(Xi*'|G)
  ## .... double check
  ## tmp.fun<-function(i, mat1=uRD.mtx, mat2=uRD.mtx){mat1[i,]%*%t(mat2[i,])};
  ## .....Daa should = matrix(apply(sapply(1:nhap, tmp.fun, mat1=uRD.mtx, mat2=uRD.mtx), 1,sum), 6,6)
  ##-----
  ## Daz
  ##-----
  Daz.1  <-t(res * XRD.post) %*% xx.adj  ## sum_i res * E(Xi*|G) * t(xx.adj)   ### modified on May 18 for allowing xx.adj
  Daz    <-Daz.1
  ##-----
  ## Dzz 
  ##-----
  Dzz.11 <-t(res * xx.adj) %*% xx.adj  ## sum_i res * xx.adj * t(xx.adj)   ### modified on May 18, 2005 for allowing xx.adj
  Dzz    <-Dzz.11
  ##-----
  ## Iaz
  ##-----
  Iaz.1   <- t(XRD.post) %*% (xx.adj * V) ## it is equal to apply(v * XRD.post,2, sum) ## modified on May 18, 2005 for xx.adj
  Iaz     <- Iaz.1
  ## Izz
  Izz.11 <- t(xx.adj * V) %*% xx.adj    ### modified on May 18 for allowing xx.adj
  Izz    <- Izz.11
  return(list("Daa"=Daa, "Daz"=Daz, "Dzz"=Dzz, "Iaz"=Iaz, "Izz"=Izz))
}





### the function below is for data analysis (allowing for missing data). Will calculate eme and determine cluster base within the fn

haplo.score.RD.unphased.fun<-function(y=Y.vec, geno=Geno, trait.type = "binomial",
                                      miss.val=NA, locus.label=NA, offset = NA, x.adj= NA, skip.haplo = 1e-7){
      ### !!!!! warning: this program can do binary and continuous Y but does not allow ordinal Y yet...
      trait.int <- charmatch(trait.type, c("gaussian", "binomial",  "poisson", "ordinal"))
      if (is.na(trait.int))        stop("Invalid trait type")
      if (trait.int == 0)          stop("Ambiguous trait type")
      if (length(y) != nrow(geno)) stop("Dims of y and geno are not compatible")
      n.loci <- ncol(geno)/2
      if (n.loci != (floor(ncol(geno)/2)))  stop("Odd number of cols of geno")
      adjusted <- TRUE
      if (all(is.na(x.adj)))       adjusted <- FALSE
      if (adjusted)
        {
          x.adj <- as.matrix(x.adj)
          if (nrow(x.adj) != length(y)) stop("Dims of y and x.adj are not compatible")
        }
      miss <- is.na(y)
      if (adjusted)       miss <- miss | apply(is.na(x.adj), 1, any)
      if (trait.int == 3)
        {
          if (all(is.na(offset))) stop("Missing offset")
          miss <- miss | is.na(offset)
        }
      y    <- as.numeric(y[!miss])
      geno <- geno[!miss, ]
      if (adjusted)       x.adj <- x.adj[!miss, , drop = FALSE]
      if (trait.int == 3) offset <- offset[!miss]
      ##-----------------------------------------------------
      ## get EME of hap freq and determine clustering base
      ##-----------------------------------------------------
      haplo <- haplo.em(geno,  locus.label=locus.label, miss.val=miss.val, control = haplo.em.control())
      H.names<-apply(haplo$haplotyp, 1, paste,  collapse="")
      Pcscn  <-haplo$hap.prob; names(Pcscn)<-H.names
      pcscn  <-rev(sort(Pcscn))
#      pos    <-getcut.fun(pcscn,nhap, plot=1)
      pos    <-getcut.fun(pcscn,nhap, plot=0)	### Modified by Wei-Chen Chen 2009-09-16
      reserv <-names(pcscn)[1:pos]
      ## modified on May 18 2005; unmark to allow for missing data

      if (!haplo$converge) stop("EM for haplo failed to converge")
      rows.rem <- haplo$rows.rem
      if (length(rows.rem) > 0)
        {
          keep <- !apply(outer((1:length(y)), rows.rem, "=="), 1, any)
          y <- y[keep]
          if (adjusted)          x.adj  <- x.adj[keep, , drop = FALSE]
          if (trait.int == 3)    offset <- offset[keep]
        }
      if (trait.int == 2)
        {
          if (!all(y == 1 | y == 0))     stop("Invalid y values")
          if (all(y == 1) | all(y == 0)) stop("No variation in y values")
        }

       ##   if (trait.int == 4) {
       ##     y     <- factor(y)
       ##     y.lev <- levels(y)
       ##     y     <- as.numeric(y)
       ##     if (max(y) < 3)             stop("Less than 3 levels for y values")
       ##   }

      ##---- prepare to get X, X*,  E(X*|G), E(X|G) ----
      n.subj     <- length(y)
      hap1       <- haplo$hap1code
      hap2       <- haplo$hap2code
      indx       <- haplo$indx.subj
      post       <- haplo$post
      nreps      <- as.vector(haplo$nreps)
      uhap       <- sort(unique(c(hap1, hap2)))
      which.haplo<- haplo$hap.prob >= skip.haplo
      p.cscn     <- Pcscn[which.haplo]
      h.names    <- names(p.cscn)
      uhap       <- uhap[which.haplo]
      ##---- obtain bigMatB for Reducing Dim ----
      BigMatB    <- final.BigMatB.matC.fun(p.cscn, reserv, n.loci, digit) ##"BB.cscn""CC.cscn""BigMatB.cscn""subPI.cscn"  
      matB       <- BigMatB$BigMatB.cscn
      CC         <- BigMatB$CC.cscn
      BB         <- BigMatB$BB.cscn
      subPI      <- BigMatB$subPI.cscn
      p.cscn.FD  <- p.cscn[match(rownames(matB), h.names)]
      p.cscn.RD  <- p.cscn.FD %*% matB
      RR         <- nrow(matB)
      RRstar     <- ncol(matB)
      ##---- get E(X*|G) =  xRD.post ----
      x          <- 0.5* (outer(hap1, uhap, "==") + outer(hap2, uhap, "==")); colnames(x)<-h.names
      x          <- x[,match(rownames(matB), h.names)]
      ## the "0.5" is because in Schaid's coding the rowsum of X vec is 2 while mine is 1
      ## xRD        <- x %*% matB[,-1]  ## modified on May 17, 2005
      xRD        <- x %*% matB
      n.xRD      <- ncol(xRD)
      xRD.post   <- matrix(rep(NA, n.subj * n.xRD), ncol = n.xRD);for(j in 1:n.xRD){ xRD.post[,j]<-tapply(xRD[,j]*post,indx,sum)}
      ##---- get E(X|G) =  x.post ----
      n.x        <- ncol(x)
      x.post     <- matrix(rep(NA, n.subj * n.x), ncol = n.x);for(j in 1:n.x){ x.post[,j]<-tapply(x[,j]*post,indx,sum)}; 
      ##-----------------------------------
      ## Note E(X*|G) should = E(X|G) %*% matB[,-1], i.e., x.post %*% matB[,-1] - xRD.post should be a 0 mat
      ##-----------------------------------
      ##     get score and Var(score)     
      ##-----------------------------------
    if (trait.int == 4) {
      stop("sorry, the program does not allow ordinal Y yet...")
       ## if (adjusted) {
       ##     library("Design")
       ##     library("Hmisc")
       ##     reg.out <- lrm(y ~ x.adj)
       ##     K <- max(y)
       ##     n.xadj <- ncol(x.adj)
       ##     alpha <- reg.out$coef[1:(K - 1)]
       ##     beta <- reg.out$coeff[K:(K - 1 + n.xadj)]
       ##     tmp <- haplo.score.podds(y, alpha, beta, x.adj, nreps, 
       ##         x.post, post, x)
       ## }
       ## if (!adjusted) {
       ##     tbl <- table(y)
       ##     s <- 1 - (cumsum(tbl) - tbl)/n.subj
       ##     alpha <- -log((1 - s[-1])/s[-1])
       ##     tmp <- haplo.score.podds(y, alpha, beta = NA, x.adj = NA, 
       ##         nreps, x.post, post, x)
       ## }
       ## u.score <- tmp$u.score
       ## v.score <- tmp$v.score
    } else  if (trait.int <= 3) {
      if (!adjusted)
        {
          mu    <- switch(trait.int, mean(y), mean(y), sum(y)/sum(offset))
          a     <- switch(trait.int, var(y), 1, 1)
          x.adj <- matrix(rep(1, n.subj), ncol = 1)
        }
      if (adjusted)
        {
          reg.out <- glm(y ~ x.adj, family = trait.type)
          x.adj   <- cbind(rep(1, n.subj), x.adj)
          mu      <- reg.out$fitted.values
          a       <- switch(trait.int, sum(reg.out$residuals^2)/reg.out$df.residual, 1, 1)
        }
      v           <- switch(trait.int, 1/a, mu * (1 - mu), mu)   ## schaid's "v"=b"/a instead of b" * a
      ##----- get score -------------       
      u.mtx     <- (y - mu) * xRD.post/a
      u.score   <- apply(u.mtx, 2, sum)  ##<-- by sum_i (yi-mu)/a * E(Xi*|G)
      ## --------------------------
      ## the result here should match with using  sum_i (yi-mu)/a * E(Xi|G) %*% matB[,-1]
      ## --------------------------
      ##  uFD.mtx  <-(y-mu)/a * x.post; uFD.score<-apply(uFD.mtx, 2, sum ;   uRD.score<- uFD.score %*% matB[,-1] 
      ## --------------------------
      ##---- get variance ----
      ## Info.mat<-get.matD.matI.fun(yy=y, mmu=mu, aa=a, XRD.post=xRD.post, XFD.post=x.post, XRD=xRD, XFD = x,
      ##                             post.prob=post, V=v, Nreps=nreps, CC.lst=CC, subPI.lst=subPI, BB.lst=BB, pFD=p.cscn.FD, xx.adj=x.adj)
      ## modified on May 19 2005. Use hybrid matEI instead of matI to allow continuous and binary Y:
      Info.mat<-get.matD.matEI.fun(yy=y, mmu=mu, aa=a, XRD.post=xRD.post, XFD.post=x.post, XRD=xRD, XFD = x,
                                  post.prob=post, V=v, Nreps=nreps, CC.lst=CC, subPI.lst=subPI, BB.lst=BB, pFD=p.cscn.FD, xx.adj=x.adj)
      ##  "Daa" "Daz" "Dzz" "Iaz" "Izz"
      if(length(Info.mat$Izz)==1){ I22i = 1/Info.mat$Izz }else{      I22i   <-Ginv(Info.mat$Izz)$Ginv }
      v.score<-(Info.mat$Daa - Info.mat$Iaz %*% I22i %*% t(Info.mat$Daz) - Info.mat$Daz %*% I22i %*% t(Info.mat$Iaz)
                + Info.mat$Iaz %*% I22i %*% Info.mat$Dzz %*% I22i %*% t(Info.mat$Iaz))
    }
      tmp   <- Ginv(v.score)
      df    <- tmp$rank
      g.inv <- tmp$Ginv
      score.global <- u.score %*% g.inv %*% u.score
      score.haplo  <- u.score/sqrt(diag(v.score))
      score.max    <- max(score.haplo^2, na.rm = TRUE)
      
      score.global.p <- 1 - pchisq(score.global, df)
      score.haplo.p  <- 1 - pchisq(score.haplo^2, 1)
      if (all(is.na(locus.label))) {   locus.label <- paste("loc-", 1:n.loci, sep = "")      }
      obj <- (list(score.global = score.global,
                   df = df,
                   score.global.p = score.global.p, 
                   score.haplo = score.haplo, 
                   score.haplo.p = score.haplo.p,
                   ## haplotype = colnames(matB)[-1],  ## modified on May 17, 2005
                   haplotype = colnames(matB),
                   matB = matB,
                   hap.prob.RD = p.cscn.RD,
                   hap.prob.FD = p.cscn.FD,
                   locus.label = locus.label
                   ))
    }




