library(rbenchmark)

# Standard way of calculating (D + UCV)^-1
classicInv <- function(A, X, sigma2){
	solve(diag(A) + (sqrt(1/sigma2)*t(X) %*% (sqrt(1/sigma2)*X)))
}


# Cholesky inverse of (D + UCV)^-1
cholInv <- function(A, X, sigma2){
	chol2inv( chol( diag(A) + (sqrt(1/sigma2)*t(X) %*% (sqrt(1/sigma2)*X))  )  )
}


# Woodbory Identity of (D + UCV)^-1
newInv <- function(A, X, sigma2){
	AInv <- diag(1/A)
	V <- sqrt(1/sigma2) * (X)
	# U = t(V)
	Imat <- diag(nrow(X))
	AInv - AInv %*% t(V) %*% solve(Imat + (V %*% AInv %*% t(V))) %*% V %*% AInv
}



# Case: n > p
set.seed(100)
n <- 1000
p <- 100
X <- matrix(rnorm(n*p), n, p)
A <- runif(p)
sigma2 <- 1.5 


system.time(classicInv(A, X, sigma2))
system.time(cholInv(A, X, sigma2))
system.time(newInv(A, X, sigma2))

benchmark(
	"classic" = {classicInv(A, X, sigma2)}, 
	"cholesky" = {cholInv(A, X, sigma2)}, 
	"woodbory" = {newInv(A, X, sigma2)}, 
	replications = 10
)




# Case: p > n 
set.seed(100)
n <- 100
p <- 1000
X <- matrix(rnorm(n*p), n, p)
A <- runif(p)
sigma2 <- 1.5 


system.time(classicInv(A, X, sigma2))
system.time(cholInv(A, X, sigma2))
system.time(newInv(A, X, sigma2))


benchmark(
	"classic" = {classicInv(A, X, sigma2)}, 
	"cholesky" = {cholInv(A, X, sigma2)}, 
	"woodbory" = {newInv(A, X, sigma2)}, 
	replications = 10
)




# Case: varying n/p ratio
p <- 1000
ratios <- c(0.01, 0.025, 0.05, 0.1, 0.25, 0.5, 1)		# --> n = (10, 25, 50, 100, 250, 500, 1000)
A <- runif(p)
sigma2 <- 1.5


for(r in ratios){

	n <- as.integer(p * r)
	print( sprintf("Ratio = %.3f", r) ) 


	set.seed(100)
	X <- matrix(rnorm(n*p), n, p)

	print(
		benchmark(
			# "classic" = {classicInv(A, X, sigma2)}, 
			"cholesky" = {cholInv(A, X, sigma2)}, 
			"woodbory" = {newInv(A, X, sigma2)}, 
			replications = 10
		)
	)
	
	print( sprintf("################################") )
	cat("\n")	

}
