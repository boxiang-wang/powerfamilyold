######################################################################################
#   All rights are reserved by the authors.
#   Authors:    Li Wang, University of Michigan (wang@umich.edu)
#               Ji Zhu, University of Michigan (jizhu@umich.edu)              
#   Date:    08/21/2004
######################################################################################
#  Parameters:
#       x: n by p matrix, where n is the data point number
#           and p is the number of predictors.
#       y: n by 1 matrix, whose element is either 1 or -1,
#           representing two classes.
#       lambda2: parameter for L2 penalize.
#       eps:    a small real number to control numberical problem
#       smallmove:  ( default value 1 ) a small positive value to
#               initialize the procedure 
#       max_steps:  the max number of iterations
#       scale:  0 or 1, where 1 means standardizing x.
#  Return:
#       beta:   s by p matrix, where s (number of steps) depends on 
#               how far the algorithm goes and p is the number of predictors
#       beta0:  s by 2 matrix (vector)
#       residuals:  length n array, where n is the number of data
#               points
#       Cvec:   length p array, which stores the last "generalized
#               correlation" for each predictor
#       loss:   length s array, which stores loss at each step
#       cor:    s by p matrix, which stores "generalized correlation" at each step.
######################################################################################
# Simple Example:
#  Training:
#       x = matrix(rnorm(80*10),nrow=80)
#	  y = sign(x[,1]+x[,2]+x[,3]+x[,4])
#	  g = svmL1L2(x,y,2)
#  Testing:
# 	  newx = matrix(rnorm(40*10),nrow=40)
#	  newy = sign(x[,1]+x[,2]+x[,3]+x[,4])
#	  pre = svmL1L2.predict(g,newx,newy)
#  Plotting:
#	  np = dim(x)
#	  n = np[1]
#	  p = np[2]
#	  horizon1 <- apply(abs(g$beta), 1, sum)
#	  matplot(horizon1, cbind(g$beta,pre$error), type="n", xlab="|beta|", ylab="beta")
#	  for(i in 1:p){
#	  	lines(horizon1, g$beta[,i], col=i+1, type="l")
#	  }
#	  lines(horizon1,pre$error,col=p+2,type="l")
######################################################################################


DrSVM_Fix2 <- function(x,y,lambda2,eps=10^(-9),smallmove=1,max_steps,scale=F) {
	np = dim(x)
	n = np[1]
	p = np[2]
	maxvars = min(p,n-1)
	lengY = length(y)
	
	if (n!=lengY){
		cat("Error:  number of rows for x and y do not match \n")
		return(NULL)
	}

	test = (abs(y) == rep(1,lengY))
	if (sum(test) != lengY){
		cat("Error:  y should only contain 1 or -1 \n")
		return(NULL)
	}
	# insert default value
	if (missing(max_steps)){
		max_steps = n*maxvars
	}

	# check invalid arguments
	if (lambda2 < 0){
		cat("Error: lambda may not be negative \n")	
		return(NULL)
	}

	if (smallmove <= 0){
        	cat("Error: smallmove must be positive \n");
        	return(NULL)
    	}

	if ((scale!=TRUE) & (scale!=FALSE)){
		cat("Error:  value for scale must be FALSE or TRUE \n");
        	return(NULL)
	}
	
	if (max_steps <=0){
		cat("Error:  value for max_steps must be greater than 0 \n")
		return(NULL)
	}
	
	indp = 1:p
	indn = 1:n
	one = rep(1,n)
	
	# standardize x
	if (scale == TRUE){
		meanx = drop(one%*%x)/n
		x = scale(x,meanx,F)	# centers x
		normx = sqrt(drop(one%*%(x^2)))
		names(normx) = NULL
		x = scale(x,F,normx)	# scales x
	}
	
	z = y*x
	
	if (sum(y==1)==sum(y==-1)){
		# balanced case
		cat("Initialize (balanced case) ...")	
		g = svmL1L2_init_bal(x,y,lambda2,eps)
		if (is.null(g)){
			cat("initial function failed\n")
			return(NULL)
		}	
		size = dim(g$beta)
		step = size[1]
		
		
		lambda1Array=g$lambda1
		
	}else{
		# unbalanced case
		cat("Initialize (unbalanced case): calling for simplex() ...")
		g = svmL1L2_init_unbal(x,y,lambda2,eps,smallmove)
		if (is.null(g)){
			cat("initial function failed\n")
			return(NULL)
		}	
		step = 1
		lambda1Array=g$lambda1

	}
	cat("done!\n")
	
	beta0 = g$beta0
	beta = g$beta
	actvar = g$actvar
	actsign = g$actsign
	correlation = g$correlation
	lambda1 = g$lambda1[step]
	
	residuals = g$residuals
	loss = g$loss
	alpha = g$alpha
	indexL = g$indexL
	indexE = g$indexE
	indexR = g$indexR

	newL = rep(FALSE,n)
	newR = rep(FALSE,n)
	newDropVar = rep(FALSE,p)

	while ((step<max_steps)&&(lambda1>0)){

		step = step + 1
		cat("step",step,":\n")
		
		if (step==2){	# unbalanced case
			beta0_val = beta0[1]
			beta_val = beta
			cor_val = correlation
		}else{	# balanced case
			beta0_val = beta0[step-1,1]
			beta_val = beta[step-1,]
			cor_val = correlation[step-1,]
		}
	
		if (sum(indexE)>0){
			len = sum(actvar) + sum(indexE)+1
			A = matrix(rep(0,len*len),nrow=len)
			b = rep(0,len)
			v = sum(actvar)
			e = sum(indexE)
			A[1:v,2:(v+1)] = lambda2 * diag(v)
			A[1:v,(v+2):(v+e+1)] = -t(z[indexE,actvar])
			A[(v+1),(v+2):(v+e+1)] = y[indexE]
			A[(v+2):(v+e+1),1] = 1
			A[(v+2):(v+e+1),2:(v+1)] = x[indexE,actvar]
			b[1:v] = actsign[actvar]			

			tmpqr = qr(A)
			if (tmpqr$rank < len){
				cat("'A is singular \n")
				step = step -1
				break
			}else{
				sol = qr.solve(tmpqr,b)
			}

			beta0_der = sol[1]
			beta_der = sol[2:(v+1)]
			alpha_der = sol[(v+2):(v+e+1)]
	
			nonezero_alphader = abs(alpha_der) > eps
			alpha_der[!nonezero_alphader] = 0	
			nonezero_betader = abs(beta_der)>eps
			beta_der[!nonezero_betader] = 0			
		}else{
			# |e| = 0
			cat("\t ------------------- |E|==0 at this step -------------------\n")
			beta0_der=0
            	beta_der = actsign[actvar]/lambda2
            	alpha_der = 0
            	nonezero_alphader = FALSE
            	nonezero_betader = abs(beta_der)>0

		}	
		# calcute delta_lam for different events
        	
		# event 1, a point reach elbow from L or R
		if (sum(actvar)>1){
			delta_fit = beta0_der + x[,actvar]%*%beta_der
		}else{
			delta_fit = beta0_der + x[,actvar]*beta_der
		}
		delta_yf = y*delta_fit
		nonezero_deltayf = abs(delta_yf) > eps
		
		if (sum(actvar)>1){
			fit = beta0_val + x[,actvar]%*%beta_val[actvar]
		}else{
			fit = beta0_val + x[,actvar]*beta_val[actvar]
		}
		temp = 1 - y*fit
			
		if (sum((!indexE)&nonezero_deltayf)>0){

			# to check if "cycling condition" exists
			if (sum(newL)>0){
				moveLeft = (delta_yf <= 0) & newL
				if (sum(moveLeft)!=sum(newL)){
					cat("\t Warning: point",indn[(!moveLeft)&newL],"just came to Left, but is coming back to Elbow!\n")
					
				}
			}			
			if (sum(newR)>0){
				moveRight = (delta_yf >= 0) & newR

				if (sum(moveRight)!=sum(newR)){
					cat("\t Warning: point",indn[(!moveRight)&newR],"just came to Right, but is coming back to Elbow!\n")
					
				}
			}

			pointOldLR = (!indexE)&(!newL)&(!newR)
			# the points on elbow do not change their residuals
			diff = temp[pointOldLR&nonezero_deltayf] / delta_yf[pointOldLR&nonezero_deltayf]
			
			if (sum(diff>0)>0){
				delta_lam1 = min(diff[diff>0])
			}else{
				delta_lam1 = Inf
			}

		}else{
			delta_lam1 = Inf
		}
		newL[] = FALSE
		newR[] = FALSE

		# event 2, a point at elbow moves to L or R (alpha moves to 1 or 0,respectively)
		elbow_alpha = alpha[indexE]
		if (sum(nonezero_alphader)>0){
			valid_alpha = elbow_alpha[nonezero_alphader]
			valid_alpha_der = alpha_der[nonezero_alphader]

			# alpha = 1 at L, alpha = 0 at R
			move_right = (valid_alpha_der < 0)
			move_left = (valid_alpha_der > 0)

			if (sum(move_right)>0){
				temp1 = valid_alpha[move_right] / (-valid_alpha_der[move_right])	# reach R
				temp1 = min(temp1)
			}else{
				temp1 = Inf
			}

			if (sum(move_left)>0){
				temp2 = (valid_alpha[move_left]-1) / (-valid_alpha_der[move_left]) # reach L
				temp2 = min(temp2)
			}else{
				temp2 = Inf
			}
			
			delta_lam2 = min(temp1,temp2)

		}else{
			delta_lam2 = Inf
		}		
		
		# event 3, a predictor reaches 0
		act_beta_val = beta_val[actvar]
		if (sum(nonezero_betader)>0){
			diff = -act_beta_val[nonezero_betader]/beta_der[nonezero_betader]
			if (sum(diff>0)>0){
				delta_lam3 = min(diff[diff>0])
			}else{
				delta_lam3 = Inf
			}

		}else{
			delta_lam3 = Inf
		}
		
		# event 4
		if (sum(actvar)<p){
			# a predictor joints {v}
			varOldInact = (!actvar)&(!newDropVar)
			dropTemp = Inf
			if (sum(newDropVar)>0){
				# some variable was just dropped, so its |correlation value|=lambda1 at this time
				# either dropTemp1 or dropTemp2 is zero
				if (sum(indexE)>1){
					dropTemp1 = (cor_val[newDropVar]+lambda1)/( t(z[indexE,newDropVar]) %*%alpha_der + 1)
					dropTemp2 = (cor_val[newDropVar]-lambda1)/(t(z[indexE,newDropVar])%*%alpha_der - 1)
				}else{
					dropTemp1 = (cor_val[newDropVar]+lambda1)/( t(z[indexE,newDropVar])*alpha_der + 1)
					dropTemp2 = (cor_val[newDropVar]-lambda1)/(t(z[indexE,newDropVar])*alpha_der - 1)

				}		
		
				if (abs(dropTemp1)<eps){
					dropTemp = dropTemp2
				}else if (abs(dropTemp2)<eps){		
					dropTemp = dropTemp1
				}else{
					cat("Error: no dropTemp is near 0 \n")
					dropTemp = Inf
				}
				if (dropTemp<0){
					dropTemp = Inf
				}

			}
			if (sum(indexE)>0){
				if (sum(varOldInact)>0){
					if (sum(indexE)>1){
						temp1 = (cor_val[varOldInact]+lambda1)/( t(z[indexE,varOldInact]) %*%alpha_der+1 )
						temp2 = (cor_val[varOldInact]-lambda1)/(t(z[indexE,varOldInact])%*%alpha_der-1)
					}else{
						temp1 = (cor_val[varOldInact]+lambda1)/( t(z[indexE,varOldInact]) *alpha_der+1 )
						temp2 = (cor_val[varOldInact]-lambda1)/(t(z[indexE,varOldInact])*alpha_der-1)

					}
					if (sum(temp1>0)>0){
						temp1 = min(temp1[temp1>0])
					}else{
						temp1 = Inf
					}

					
					if (sum(temp2>0)>0){
						temp2 = min(temp2[temp2>0])
					}else{
						temp2 = Inf
					}

					delta_lam4 = min(c(temp1,temp2,dropTemp))
					if (delta_lam4 > lambda1){
						delta_lam4 = Inf
					}	
				}else{

					delta_lam4 = dropTemp
				}
			}else{
				delta_lam4 = min(lambda1 - abs(cor_val[varOldInact]))
			}
			

		}else{
			# all predictors are included    
            	delta_lam4 = lambda1;   # the distance to have all correlations reach zero
		}
		
		newDropVar[] = FALSE

		seq = c(delta_lam1,delta_lam2,delta_lam3,delta_lam4)

		delta_lam = min(seq)

		if (sum(seq==delta_lam)>1){
			cat("two events tie abnormally. \n")
			step = step - 1
			break
		}

		if (delta_lam>lambda1){
			cat("\t Stop: no further events can happen, lambda1 will decrease to 0\n")
			delta_lam = lambda1
		}

		beta0_val = beta0_val + beta0_der*delta_lam
		beta_val[actvar]=beta_val[actvar]+beta_der*delta_lam
		
		if (sum(indexE)>0){
			alpha[indexE] = alpha[indexE]+alpha_der*delta_lam
		}

		if (sum(actvar)>1){
			fit = beta0_val+x[,actvar]%*%beta_val[actvar]
		}else{
			fit = beta0_val+x[,actvar]*beta_val[actvar]
		}
		
		temp = 1 - y*fit
		residuals = temp
		residuals[residuals<0]=0
		lambda1 = lambda1 - delta_lam
		cor_val = lambda2*beta_val - alpha%*%z
		
		cat("\t Generalized Corr -",delta_lam,"to",lambda1)
		
		# update actvar, actsign, indexL,indexE and indexR		

		if (delta_lam==delta_lam1){
			# event 1, a point reach elbow from L or R - update indexL,indexE and indexR
			
			index_fromL = (abs(temp)<eps) & indexL
			index_fromR = (abs(temp)<eps) & indexR			

			if (sum(index_fromL)>0){
				cat("\t point",indn[index_fromL],"reached Elbow from Left\n")
				indexL[index_fromL] = FALSE
				indexE[index_fromL] = TRUE
				residuals[indexE]=0
			}

			if (sum(index_fromR)>0){
				cat("\t point",indn[index_fromR],"reached Elbow from Right\n")
				indexR[index_fromR] = FALSE
				indexE[index_fromR] = TRUE
				residuals[indexE]=0
			}

			if (sum(index_fromL|index_fromR)==0){
				cat("\t numerical error 1\n")
			}
			
			if (sum(actvar)<(sum(indexE)-1)){
				cat("Warning: |v| < |e| - 1 \n")
			}

		}else if (delta_lam==delta_lam2){
			# a point at elbow moves to L or R (residuals still 0, alpha moves to 1 or 0,respectively)
						
			valid_alpha_der = alpha_der[nonezero_alphader]
			move_right = valid_alpha_der < 0
			move_left = valid_alpha_der > 0
			if (sum(nonezero_alphader)==0){
				cat("Error: alpha_der == 0\n")
			}

			tempIndexE = indexE
			tempIndexE[tempIndexE] = nonezero_alphader
			index_toL = (abs(alpha-1)<eps)&tempIndexE
			index_toL[tempIndexE] = index_toL[tempIndexE] & move_left

			index_toR = (abs(alpha)<eps)&tempIndexE
			index_toR[tempIndexE] = index_toR[tempIndexE] & move_right

			if (sum(index_toL)>0){
				cat("\t point",indn[index_toL],"left Elbow to the Left\n")
				indexL[index_toL] = TRUE
				indexE[index_toL] = FALSE
				newL[index_toL] = TRUE
				alpha[index_toL] = 1
			}
			
			if (sum(index_toR)>0){
				cat("\t point",indn[index_toR],"left Elbow to the Right\n")
				indexR[index_toR] = TRUE
				indexE[index_toR] = FALSE
				newR[index_toR] = TRUE
				alpha[index_toR] = 0
			}
			if (sum(index_toL|index_toR)==0){
				cat("\t numerical error 2\n")
				break
			}
			
		}else if (delta_lam==delta_lam3){
			# event 3, a predictor reaches 0
			index_toZero = (abs(beta_val)<eps)&actvar
			
			if (sum(index_toZero)>0){
				cat("\t predictor",indp[index_toZero],"left {V} \n")
				actvar[index_toZero] = FALSE
				actsign[index_toZero] = 0
				beta_val[index_toZero] = 0
				newDropVar[index_toZero] = TRUE

			}else{
				cat("\t numerical error 3\n")
				break
			}
			
			if (sum(actvar)<(sum(indexE)-1)){
				cat("Warning: |v| < |e| - 1\n")
				break
			}

		}else if (delta_lam==delta_lam4){
			if (sum(actvar)<p){
				# a predictor joins {v}
				index_toV = (abs(abs(cor_val)-max(abs(cor_val)))<eps)&(!actvar);
				if (sum(index_toV)>0){
					cat("\t predictor",indp[index_toV],"joined {V}\n")
					actvar[index_toV] = TRUE
					actsign[index_toV] = -sign(cor_val[index_toV])
				}else{
					cat("\t numerical error 4\n")
					break
				}
			} # otherwise all predictors are in {v} and corr will reach 0 in the next step 
			  # it will exit the while loop next step 

		}else{
			cat("error: delta_lam has abnormal value \n")
			break
		}
		
		beta0 = rbind(beta0,c(0,0))
		beta0[step,] = beta0_val

		beta = rbind(beta,rep(0,p))
		beta[step,] = beta_val

		correlation = rbind(correlation,rep(0,p))
		correlation[step,] = cor_val

		loss = c(loss,0)
		loss[step] = sum(residuals)
	
		lambda1Array[step]=lambda1

		max_cor = max(abs(cor_val))

		if (loss[step]<eps){
			cat("\t Stop: loss is 0 now!\n")
			break
		}else if (max_cor<eps){
			cat("\t Stop: all correlations reach 0 now!\n")
			break
		}

	}	# while
	if (step == max_steps){
		cat("\t Stop: max_steps is reached now!\n")
	}
	Cvec = cor_val
	cor = correlation

	return(list(beta=beta,beta0=beta0,residuals=residuals,
	Cvec=Cvec,loss=loss,cor=cor,lambda1=lambda1Array))
}

svmL1L2_init_bal <- function(x, y, lambda2, eps=10^(-9)) {
### y \in {1,-1}
	np = dim(x)
	n = np[1]
	p = np[2]
	if (sum(y==1)!=sum(y==-1)){
		cat("Error: y is not balanced \n")
		return(NULL)
	}
	step = 1
	beta0 = c(-1,1)	#initial range of beta0
	beta = rep(0,p)	#initial value of beta

	#Now alpha_i = 1 for all points
	cor = y %*% x
	jstar = (abs(cor)==max(abs(cor)))
	
	if (sum(jstar)>1){
		cat("Error: more than one predictor enters abnormally at 1\n")
		return(NULL)
	}
	
	sign_beta = sign(cor[jstar])	# it's a scalor 
	lam_init = max(abs(cor)) # when lambda1 decrease from Inf to lam_init, the first predictor joins
	
	lambda1 = lam_init

	correlation = lambda2*beta - cor

	beta_der = sign_beta / lambda2	#  derivative of beta (d_beta/d_lambda1) - scalor

	cx = beta_der * x[,jstar]
	i_plus = (cx == max(cx[y==1]))
	i_neg = (cx == min(cx[y==-1]))
	diff = cx[i_plus] - cx[i_neg]
	if (diff <=0){
		cat("Abnormal data: cx(i_plus) < cx(i_neg) \n")
		return(NULL)
	}
	
	delta_lambda = 2/diff
	temp_cor = cor
	temp_cor[jstar] = 0
	jstar2 = (abs(temp_cor)==max(abs(temp_cor))) # second largest correlation
	if (sum(jstar2)>1){
		cat("more than one predictor enters abnormally at 2\n")
		return(NULL)
	}
	lam_sec = max(abs(temp_cor))
	delta_lambda2 = lam_init - lam_sec	# distance needed to have another predictor join {v}
	step = step + 1

	actvar = jstar
	actder = rep(0,p)
	actder[jstar] = beta_der
	actsign = rep(0,p)
	actsign[jstar] = sign_beta
	
	uplow_bound = y
	while (delta_lambda2 < delta_lambda){ # new predictor joint {v} before beta0 is fixed
		actvar = actvar | jstar2
		actsign[jstar2] = sign(cor[jstar2])
		lam_init = lam_sec
		
		lambda1[step] = lam_sec

		beta0 <- rbind(beta0,c(0,0))

		uplow_bound[y==1] = uplow_bound[y==1] - delta_lambda2*cx[y==1]
		uplow_bound[y==-1] = uplow_bound[y==-1] - delta_lambda2*cx[y==-1]
		
		beta0[step,1] = max(uplow_bound[y==-1]) # updated largest lower bound for beta0
		beta0[step,2] = min(uplow_bound[y==1]) # updated smallest upper bound for beta0
			
		beta = rbind(beta,rep(0,p))
		beta[step,] = beta[step-1,] + actder * delta_lambda2 # updated beta
		correlation = rbind(correlation,rep(0,p))
		correlation[step,] = lambda2*beta[step,] - cor
		
		actder = actsign / lambda2
		# compute delta_lambda, the distance needed to have beta0 fixed
		if (sum(actvar)>1){
			cx = x[,actvar]%*%actder[actvar]
		}else{
			cx = x[,actvar]*actder[actvar]
		}
		res_plus = uplow_bound[y==1]
		res_neg = uplow_bound[y==-1]
		cx_plus = cx[y==1]
		cx_neg = cx[y==-1]
		one_matrix = matrix(rep(1,n*n/4),nrow=n/2)		
		
		delta_matrix = res_plus*one_matrix - t(res_neg*one_matrix)
		delta_matrix = delta_matrix / (cx_plus*one_matrix - t(cx_neg*one_matrix))
		delta_matrix[delta_matrix <0]=Inf
		delta_lambda = min(delta_matrix)
		if (delta_lambda == 0){
			cat("Error: delta_lambda = 0 \n")
			return(NULL)
		}
		
		# compute delta_lambda2, the distance needed to have anther
        	# predictor join {v}
        	temp_cor=cor
        	temp_cor[actvar] = 0
		jstar2 = (abs(temp_cor)==max(abs(temp_cor)))	# find the next largest correlation
		if (sum(jstar>1)){
			cat("more than one predictor enters abnormally at",step+1,"\n")
			return(NULL)
		}
		lam_sec = max(abs(temp_cor))
		delta_lambda2 = lam_init - lam_sec
		
		step = step + 1
	}	# while
	if (delta_lambda2 == delta_lambda){
		cat("delta_lambda ties delta_lambda2 abnormally \n");
		return(NULL)
	}
	uplow_bound[y==1] = uplow_bound[y==1] - delta_lambda*cx[y==1]
	uplow_bound[y==-1] = uplow_bound[y==-1] - delta_lambda*cx[y==-1]
	
	beta0 = rbind(beta0,rep(0,2))
	beta0[step,1] = max(uplow_bound[y==-1]) # updated largest lower bound for beta0
	beta0[step,2] = min(uplow_bound[y==1])  # updated smallest upper bound for beta0
	if (abs(beta0[step,1]-beta0[step,2]) > eps){
		cat("numerical error: beta0 is not fixed \n")
	}
	
	beta = rbind(beta,rep(0,p))
	beta[step,] = beta[step-1,] + actder*delta_lambda
	actsign = sign(beta[step,])

	correlation = rbind(correlation,rep(0,p))
	correlation[step,] = lambda2*beta[step,] - cor
	lambda1[step] = max(abs(correlation[step,]))
	
	if (sum(actvar)>1){
		fit = beta0[step,1] + x[,actvar] %*% beta[step,actvar]
	}else if (sum(actvar)==1){
		fit = beta0[step,1] + x[,actvar] * beta[step,actvar]
	}else{
		cat("error: no actvar \n")
		return(NULL)
	}
	residuals = 1- y*fit
	residuals[abs(residuals)<eps] = 0
	residuals[residuals<0] = 0
	loss = rep(0,step) + NaN
	loss[step] = sum(residuals)

	alpha = rep(1,n)
	indexL = residuals > eps
	if (length(indexL) != length(y)){
		cat("initial error: all points should be on L \n")
		return(NULL)
	}
	indexE = !indexL
	indexR = indexL & indexE
	
	return(list(beta0=beta0,beta=beta,actvar=actvar,actsign=actsign,
	correlation=correlation,lambda1=lambda1,residuals=residuals,loss=loss,
	alpha=alpha,indexL=indexL,indexE=indexE,indexR=indexR))
}

svmL1L2_init_unbal <- function(x, y, lambda2, eps=10^(-9), smallmove = 1) {
	np = dim(x)
	n = np[1]
	p = np[2]	
	
	indn = 1:n
	indp = 1:p

	require(boot)
	
	if (sum(y==1) > sum(y==-1)){
		beta0 = 1
		index_bigclass = (y==1)
		index_smallclass = (y==-1)
	}else if (sum(y==1)<sum(y==-1)){
		beta0 = -1
		index_bigclass = (y==-1)
		index_smallclass = (y==1)
	}else{
		# sum(y==1) == sum(y==-1)
		cat("Error: y is balanced \n")
		return(NULL)
	}

	xbig = x[index_bigclass,]
	xsmall = x[index_smallclass,]
	ybig = y[index_bigclass]
	ysmall = y[index_smallclass]
	nbig = length(ybig)
	nsmall = length(ysmall)

	one = rep(1,nsmall)
	tmp = one %*% xsmall
	tmp = ysmall[1]*t(tmp)
	f = c(rep(1,nbig),-ysmall[1]*nsmall,-tmp,ysmall[1]*nsmall,tmp)
	s = smallmove
	A1 = matrix(c(rep(0,nbig),0,rep(1,p),0,rep(1,p)),nrow=1)
	b1 = s
	A2 = cbind(diag(nbig),ybig,ybig*xbig,-ybig,-ybig*xbig)
	b2 = rep(1,nbig)
	g = simplex(f,A1=A1,b1=b1,A2=A2,b2=b2,maxi = FALSE, n.iter= n*p)
	if (g$solved == -1){
		cat("simplex() can not converge to a solution \n")
		return(NULL)
	}else if (g$solved == 0){
		cat("Warning: simplex() reaches the maximum iterations \n")
		return(NULL)
	}

	tmp = g$soln[-(1:nbig)]
	beta0 = tmp[1] - tmp[p+2]
	beta = tmp[2:(p+1)] - tmp[(p+3):(2*p+2)]
	
	actvar = abs(beta) > eps
	actvar = t(actvar)
	beta[abs(beta) <=eps] = 0
	actsign = sign(beta)
	if (p>1){
		res = 1 - y*(beta0 + x %*% beta)	
	}else if (p==1){
		res = 1 - y*(beta0 + x * beta)
	}else{
		cat("Error: p = 0\n")
		return(NULL)
	}
	indexE = (abs(res)<=eps)&index_bigclass
	indexL = index_smallclass | ((res>eps)&index_bigclass)
	indexR = index_bigclass & (res<(-eps))
	sign_beta = sign(beta)

	if (sum(indexE)!= sum(actvar)){
		cat("Elbow points =", sum(indexE),", Num of predictors =",sum(actvar),"\n")
		if ((sum(indexE)==0) && (sum(actvar)==1)){
			minVal = min(abs(res[index_bigclass]))
			nearestPnt = (abs(res) == minVal)&index_bigclass
			cat("Warning: there is no unique solution for beta0, we force the nearest point to stay on Elbow. \n")
			
			if (sum(nearestPnt)>1){	# two points have the same abs residuals
				cat("Warning: Two points have the same absolute residual value\n")
				temp = res + runif(1)*eps # shift by a random value
				minVal = min(abs(temp[index_bigclass]))
				nearestPnt = (abs(temp)==minVal)&index_bigclass
				if (sum(nearestPnt) > 1){
					cat("Warning: Two points have exactly the same residual value\n")
					return(NULL)
				}
			}			
			indexE[nearestPnt] = 1
			indexL[nearestPnt] = 0
			indexR[nearestPnt] = 0

		}else{
			cat("Warning: simplex() does not converge to a correct minimum, try other smallmove value \n")
			return(NULL)
		}
	}
	
	beta = rep(0,p)	# initial beta, which will be returned
	beta0 = rep(ybig[1],2)  # initial beta0, which will be returned

	alpha = rep(0,n)
	alpha[indexL] = 1
	A = matrix(rep(0,(sum(actvar)+1)*(sum(indexE)+1)),nrow = sum(actvar)+1)	
	b = rep(1,sum(actvar)+1)
	b[1:sum(actvar)] = t(x[,actvar])%*% (alpha * y)
	b[sum(actvar)+1] = -sum(alpha * y)
	if (sum(indexE)>1){
		A[1:sum(actvar),1:sum(indexE)] = -t(x[indexE,actvar])%*% diag(y[indexE])
    	}else{
		A[1:sum(actvar),1:sum(indexE)] = -t(x[indexE,actvar]) * (y[indexE])
	}
	A[1:sum(actvar),(sum(indexE)+1)] = sign_beta[actvar]
    	A[(sum(actvar)+1),1:sum(indexE)] = t(y[indexE])

	tmpqr = qr(A)
	if (tmpqr$rank < length(b)){
		cat("'A is singular, alphas can not be calculated \n")
		return(NULL)
	}else{
		sol = qr.solve(tmpqr,b)
	}

	alpha[indexE] = sol[1:sum(indexE)]
	lam_init = sol[sum(indexE)+1]	# lam_init is the critical lambda1 value

	if (abs(sum(alpha * y)) > eps){
		cat("error: sum(alpha_i * y_i) is not 0\n")
		return(NULL)
	}
	if (sum(alpha<0)>0){
		if (sum(0 - alpha[alpha<0]) > eps){
			cat("error: negative alpha exists\n")	
			return(NULL)
		}else{
			alpha[alpha<0]=0
		}
	}
	if (sum(alpha>1)>0){
		if (sum(alpha[alpha>1]-1)>eps){
			cat("error: some alpha has value > 1\n")
			return(NULL)
		}else{
			alpha[alpha>1] = 1
		}
	}
	if (lam_init <= 0){
		cat("error: lambda has negative initial value \n")
		return(NULL)
	}
	correlation = - (alpha*y) %*% x
	lambda1 = max(abs(correlation))
	residuals = 1 - y*beta0[1]
	loss = sum(residuals)
	
	return(list(beta0=beta0,beta=beta,actvar=actvar,actsign=actsign,
	correlation=correlation,lambda1=lambda1,residuals=residuals,loss=loss,
	alpha=alpha,indexL=indexL,indexE=indexE,indexR=indexR))

}	


DrSVM_Fix2.predict <- function(obj, x, y=NULL) {
### obj: output of svmL1L2
### y \in {1,-1}
  N <- length(y)
  f <- obj$beta %*% t(x) + obj$beta0[,1]	# 1 and 2 produce the same fitting value
  predict <- sign(f)
  if (is.null(y))
    return(list(f=f, error=NULL))
  error <- array(apply(apply(predict, 1, FUN="!=", y), 2, sum)/N)
  return(list(f=f, error=error))
}
                            
