attach(faithful)
y = waiting
EM_gaussian = function(x, theta){
  # E-step:
  bstar = theta[1]*dnorm(x, theta[2], theta[3])/(theta[1]*dnorm(x, theta[2], theta[3])+
                                                   (1-theta[1])*dnorm(x, theta[4], theta[5]))
  
  # M-step:
  n = length(x)
  theta_M = c()
  theta_M[1] = mean(bstar)
  theta_M[2] = sum(bstar*x)/sum(bstar)
  theta_M[3] = sum(bstar*(x-theta_M[2])^2)/sum(bstar)
  theta_M[4] = sum((1-bstar)*x)/sum(1-bstar)
  theta_M[5] = sum((1-bstar)*(x-theta_M[4])^2)/sum(1-bstar)
  
  return(c(theta_M[1],theta_M[2],sqrt(theta_M[3]),theta_M[4],sqrt(theta_M[5])))
}

# Running the EM algorithm

theta = c()
bstar = NULL
theta0 = c(0.4, 52, 5, 80, 5)
dis = 1
i = 1
a0 = NA; a1 = NA
while (dis > 0.001 || i <= 200){
  
  theta = EM_gaussian(y,theta0)
  dis = max(abs(theta0 - theta))
  if(i == 199){a0 = theta - theta0}
  if(i == 200){a1 = theta - theta0}
  theta0 = theta
  i = i+1
  print(i)
}
theta0
a0
a1
a1/a0



