



set.seed(1234567)
mX1 <- rnorm(100,mean = 0, sd = 1)
mX2 <- rnorm(100,mean = 1, sd = 2)


# Parameters 
a = 0.9
b1 = 0.5
b2 = -0.2
c1 = 0.01
c2 = 0.04

y = a + (b1 * mX1) + (b2 * mX2)

