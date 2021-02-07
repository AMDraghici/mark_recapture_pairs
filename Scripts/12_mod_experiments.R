library(rjags)

model_file <- "
  model{  
    
    y[1] ~ dbern(1)
    for(i in 2:newIDn){
          y[i] ~ dbern(p[chosen_index[i]])
    }
  
    for(i in 1:newIDn){
      chosen_index[i] ~ dcat(d[1:oldIDn]) 
    }
    
   for(j in 2:oldIDn){
     is_chosen_sum[j] <- sum(chosen_index[1:newIDn]==j)
     zeroes_vec[j] ~ dinterval(is_chosen_sum[j], 1)
   }
    
    
    for(i in 1:oldIDn){
      p[i] ~ dunif(0, 1)
    }
    


    for(i in 1:oldIDn){
      d[i] <- exp(eta[i])/exp.norm[oldIDn]
    }
    
    
    exp.norm[1] <- exp(eta[1])
    for(i in 2:oldIDn){
      exp.norm[i] <- exp(eta[i]) + exp.norm[i-1]
    }
    
  
  eta[1:oldIDn] <-  beta0 + beta1*x[1:oldIDn]
    
   #prior for probs
    beta0 ~ dnorm(0,1)
    beta1 ~ dnorm(0,1)
  }
"
tmpf=tempfile()
tmps=file(tmpf,"w")
cat(model_file,file=tmps)
close(tmps)

newIDn <- 100
oldIDn <- 100

## Generate values
jags.data <- list(newIDn = newIDn,
                  oldIDn = oldIDn,
                  chosen_index = sample(oldIDn,replace=FALSE),
                  x = runif(oldIDn, -5, 5),
                  zeroes_vec = rep(0,oldIDn),
                  y = rbinom(newIDn,1,0.5))

jags.data$chosen_index[c(3,5,7)] <- NA
mymodel <- jags.model(tmpf,jags.data)
iter <- 100
mysamples <- coda.samples(mymodel,c("p","beta0","beta1"),iter)

hist(mysamples[[1]][,5])
hist(mysamples[[1]][,7])
hist(mysamples[[1]][,9])
#summary(mysamples)

#ggs_density(ggs(mysamples))
