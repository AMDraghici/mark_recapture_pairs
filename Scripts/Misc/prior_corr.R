p1 <- rbeta(1e6, 3,3)
# odds_p1 <- p1/(1-p1)
# 
# beta <- rgamma(1e6,3,1)
# alpha <- beta * odds_p1 #rgamma(1e6,(beta * odds_p1), 1)
v<- rgamma(1e6,10,1)
alpha <- p1 * v
beta <- (1-p1) * v


plot(density(beta))
plot(density(alpha))

p2 <- rbeta(1e6,alpha,beta)

cor(p1,p2)
lm(p2 ~ p1 - 1)

data.frame(p1 = p1, p2 = p2) %>% 
  ggplot(aes(x = p1, y = p2)) +
  stat_ellipse(level = 0.95) +
  stat_ellipse(level = 0.99) +
  stat_ellipse(level = 0.50) +
  geom_smooth(method = "lm")


data.frame(p1=p1,p2=p2) %>% 
  filter(p1 < 0.91 & p1 > 0.89) %>% 
  pull(p2) %>% 
  density() %>% 
  plot()

# plot(density(p1))
# lines(density(p2), col = "red")

data.frame(p1=p1,p2=p2) %>% 
  ggplot() +
  geom_density(aes(x = p1)) +
  geom_density(aes(x = p2), col = "red") +
  geom_density(aes(x = runif(1e6,0,1)), col = "blue")

