### Examples to use EM functions.
data.path <- paste(.libPaths()[1], "/phyclust/data/pony524.phy", sep = "")
my.seq <- read.phylip(data.path)
X <- my.seq$org
my.seq <- read.phylip(data.path)
X <- my.seq$org

### Directly use phyclust().
(ret <- phyclust(X, 2))

### One EM step.
(ret.em <- phyclust.em.step(X, ret))

### One E- and M- step.
ret.e <- phyclust.e.step(X, ret)
(ret.m <- phyclust.m.step(X, ret))
(ret.e.m <- phyclust.m.step(X, K = ret$K, Z.normalized = ret.e,
                            substitution.model = ret$substitution.model,
                            identifier = ret$QA$identifier,
                            code.type = ret$code.type))

### Log-likelihood.
(ret.logL <- phyclust.logL(X, ret))
