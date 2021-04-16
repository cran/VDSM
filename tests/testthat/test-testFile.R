library(VDSM)

data(exampleX)
X=exampleX
data(examplef)
f=examplef
p=8
Anchor.estimate=c(3,2.5,2,1.5,1,0,0,0)

test_that("Gplot works", {
  G_example1 = Gplot(X,f,p)
})
