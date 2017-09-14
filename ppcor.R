#Relevent code for scTarNet taken from the ppcor package 

# partial correlation
pcor <- function(x, method = c("pearson", "kendall", "spearman"))
{
  # correlation method
  method <- match.arg(method)
  
  # check the data
  if (is.data.frame(x)) 
    x <- as.matrix(x)
  if (!is.matrix(x)) 
    stop("supply a matrix-like 'x'")
  if (!(is.numeric(x) || is.logical(x))) 
    stop("'x' must be numeric")
  stopifnot(is.atomic(x))
  
  # sample number
  n <- dim(x)[1]
  
  # given variables' number
  gp <- dim(x)[2]-2
  
  # covariance matrix
  cvx <- cov(x,method=method)
  
  # inverse covariance matrix
  if(det(cvx) < .Machine$double.eps){
    warning("The inverse of variance-covariance matrix is calculated using Moore-Penrose generalized matrix inverse due to its determinant of zero.")
    library('MASS')
    icvx <- ginv(cvx)
  }else
    icvx <- solve(cvx)
  
  # partial correlation
  pcor <- -cov2cor(icvx)
  diag(pcor) <- 1
  
  # p-value
  if(method == "kendall"){
    statistic <- pcor/sqrt(2*(2*(n-gp)+5)/(9*(n-gp)*(n-1-gp)))
    p.value <- 2*pnorm(-abs(statistic))
    
  }else{
    statistic <- pcor*sqrt((n-2-gp)/(1-pcor^2))
    p.value <- 2*pt(-abs(statistic),(n-2-gp))
    #p.value <- 2*pnorm(-abs(statistic))
  }
  
  diag(statistic) <- 0
  diag(p.value) <- 0
  
  list(estimate=pcor,p.value=p.value,statistic=statistic,n=n,gp=gp,method=method)
}

# pairwise partial correlation
pcor.test <- function(x,y,z,method=c("pearson", "kendall", "spearman"))
{
  # The partial correlation coefficient between x and y given z
  #
  # pcor.test is free and comes with ABSOLUTELY NO WARRANTY.
  #
  # x and y should be vectors
  #
  # z can be either a vector or a matrix
  
  # correlation method
  method <- match.arg(method)
  
  x <- c(x)
  y <- c(y)
  z <- as.data.frame(z)
  
  # merge into a matrix
  xyz <- data.frame(x,y,z)
  
  # partial correlation
  pcor = pcor(xyz,method=method)
  
  data.frame(estimate=pcor$est[1,2],p.value=pcor$p.value[1,2],statistic=pcor$statistic[1,2],n=pcor$n,gp=pcor$gp,Method=method)
}	

