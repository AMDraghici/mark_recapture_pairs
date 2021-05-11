set.seed(1)
x <- list()
y <- list()
for(i in 1:3){
  x[[i]] <- runif(3,0,1)
  y[[i]] <- runif(3,0,1)
}

for(i in 1:3){
  x[[i]] <- x[[i]]/sum(x[[i]])
  y[[i]] <- y[[i]]/sum(y[[i]])
}

x <- do.call(rbind, x)
y <- do.call(rbind, y)

rowSums(x)
colSums(t(y))

# 1. Show how order matters (what to do)
(1-x[1,1])*(x[1,2]+(1-x[1,2])*x[1,3]+(1-x[1,2])*(1-x[1,3]))*x[2,1]
(1-x[1,1])*x[2,1]
(1-x[1,1])*(1-x[2,1])*x[3,1]

# 2. Show prob of not picking partner (solution is mostly figured out)
(1-x[1,1])*(1-x[1,2])*(1-x[1,3])
(1-y[1,1])*(1-y[1,2])*(1-y[1,3])

# Solution: condition across rows

# Prob of partner 1 (Prob (1 and 1 form a pair))

# Prob of partner 2 
##(prob[12|1!=1] = prob[12 n 1!=1]/prob[1!=1] = prob[1!=1|12]*prob[12]/prob[1!=1])
### prob[1!=1|12] = 1 
x[1,2]/(1-x[1,1])

# Prob of partner 3
##(prob[13|1!=1,1!=2] = prob[13 n 1!=1 n 1!=2]/prob[1!=1 n 1!=2] = prob[1!=1 n 1!=2|13]*prob[13]/(1 - p(1=1) - p(1=2)))
###prob[1!=1 n 1!=2|13] = 1
####prob[1!=1 n 1!=2] = P[!(1=1 U 1=2)] = 1 - P[(1=1 U 1=2)] = 1 - p(1=1) - p(1=2) (since p(1=1 n 1=2) = 0)
x[1,3]/(1-x[1,1]-x[1,2])

# theres also a more convoluted way to find the denominator using conditional probs

# if we fix problem two (basically already resolved) then the problem in part 1 becomes worse 

x <- matrix(c(1/4,3/4,1/4,3/4), byrow = T, nrow = 2, ncol =2)

# We go from 
(1-x[1,1])*x[1,2]

# to 
x[1,2]/(1-x[1,1])
(1-x[1,1])
# basically a 75% chance for pair 1 = 2 to happen
# 25% chance for pair 2-2 to happen
# actually the distribution of x[2,] is irrelevant 
