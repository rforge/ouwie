#OUwie Master Controller

#written by Jeremy M. Beaulieu

#Fits the Ornstein-Uhlenbeck model of continuous characters evolving under discrete selective 
#regimes. The input is a tree of class "phylo" that has the regimes as internal node labels 
#and a trait file. The trait file must be in the following order: Species names, Regime, and 
#continuous trait. Different models can be specified -- Brownian motion (BM), multiple rate BM (BMS)
#global OU (OU1), multiple regime OU (OUM), multiple sigmas (OUMV), multiple alphas (OUMA), 
#and the multiple alphas and sigmas (OUMVA). 

OUwie<-function(phy,data, model=c("BM1","BMS","OU1","OUM","OUMV","OUMA","OUMVA"), simmap.tree=FALSE, root.station=TRUE, ip=1, lb=0.000001, ub=1000, clade=NULL, quiet=FALSE){
	
	#Makes sure the data is in the same order as the tip labels
	data<-data.frame(data[,2], data[,3], row.names=data[,1])
	data<-data[phy$tip.label,]
	#Values to be used throughout
	n=max(phy$edge[,1])
	ntips=length(phy$tip.label)
	ip=ip
	#Will label the clade of interest if the user so chooses:
	if(is.null(clade)){
		phy=phy
	}
	if(!is.null(clade) & simmap.tree==FALSE){
		node<-mrca(phy)[clade[1],clade[2]]
		int<-c(node,Descendants(phy,node, "all"))
		tips<-int[int<ntips]
		data[,1]<-1
		data[tips,1]<-2
		int<-int[int>ntips]
		phy$node.label<-rep(1,phy$Nnode)
		pp<-int-length(phy$tip.label)
		phy$node.label[pp]<-2
	}
	if (is.character(model)) {
		if (model == "BM1"| model == "OU1"){
			simmap.tree=FALSE
		}
	}
	if(simmap.tree==TRUE){
		k<-length(colnames(phy$mapped.edge))
		tot.states<-factor(colnames(phy$mapped.edge))
		tip.states<-factor(data[,1])
		data[,1]<-as.numeric(tip.states)
		
		#A boolean for whether the root theta should be estimated -- default is that it should be.
		root.station=root.station
		#Obtains the state at the root
		root.state=which(colnames(phy$mapped.edge)==names(phy$maps[[1]][1]))
		##Begins the construction of the edges matrix -- similar to the ouch format##
		#Makes a vector of absolute times in proportion of the total length of the tree
		branch.lengths=rep(0,(n-1))
		branch.lengths[(ntips+1):(n-1)]=branching.times(phy)[-1]/max(branching.times(phy))
		
		#New tree matrix to be used for subsetting regimes
		edges=cbind(c(1:(n-1)),phy$edge,phy$edge.length)
		edges=edges[sort.list(edges[,3]),]
		
		edges[,4]=branch.lengths
		
		#Resort the edge matrix so that it looks like the original matrix order
		edges=edges[sort.list(edges[,1]),]
	}
	if(simmap.tree==FALSE){
		#Obtain a a list of all the regime states. This is a solution for instances when tip states and 
		#the internal nodes are not of equal length:
		tot.states<-factor(c(phy$node.label,as.character(data[,1])))
		
		k<-length(levels(tot.states))

		int.states<-factor(phy$node.label)
		phy$node.label=as.numeric(int.states)		
		tip.states<-factor(data[,1])
		data[,1]<-as.numeric(tip.states)
		
		#A boolean for whether the root theta should be estimated -- default is that it should be.
		root.station=root.station
		if (is.character(model)) {			
			if (model == "BM1"| model == "OU1"){
				##Begins the construction of the edges matrix -- similar to the ouch format##
				#Makes a vector of absolute times in proportion of the total length of the tree
				phy$node.label<-sample(c(1:k),phy$Nnode, replace=T)
				int.states=length(levels(tip.states))
				#Since we only really have one global regime, make up the internal nodes -- this could be improved
				phy$node.label<-as.numeric(int.states)
				branch.lengths=rep(0,(n-1))
				branch.lengths[(ntips+1):(n-1)]=branching.times(phy)[-1]/max(branching.times(phy))
				
				#Obtain root state -- for both models assume the root state to be 1 since no other state is used even if provided in the tree
				root.state<-1
				#New tree matrix to be used for subsetting regimes
				edges=cbind(c(1:(n-1)),phy$edge,phy$edge.length)
				edges=edges[sort.list(edges[,3]),]
				
				edges[,4]=branch.lengths
				regime <- matrix(0,nrow=length(edges[,1]),ncol=k)
				regime[,1]<-1
				regime[,2:k]<-0
				
				edges=cbind(edges,regime)
			}
			else{
				##Begins the construction of the edges matrix -- similar to the ouch format##
				#Makes a vector of absolute times in proportion of the total length of the tree
				branch.lengths=rep(0,(n-1))
				branch.lengths[(ntips+1):(n-1)]=branching.times(phy)[-1]/max(branching.times(phy))
				
				#Obtain root state and internal node labels
				root.state<-phy$node.label[1]
				int.state<-phy$node.label[-1]
				#New tree matrix to be used for subsetting regimes
				edges=cbind(c(1:(n-1)),phy$edge,phy$edge.length)
				edges=edges[sort.list(edges[,3]),]
				
				edges[,4]=branch.lengths
				
				mm<-c(data[,1],int.state)
				regime <- matrix(0,nrow=length(mm),ncol=length(unique(mm)))
				#Generates an indicator matrix from the regime vector
				for (i in 1:length(mm)) {
					regime[i,mm[i]] <- 1 
				}
				#Finishes the edges matrix
				edges=cbind(edges,regime)
			}
		}
		#Resort the edge matrix so that it looks like the original matrix order
		edges=edges[sort.list(edges[,1]),]
	}
	x<-as.matrix(data[,2])
	#Matches the model with the appropriate parameter matrix structure
	if (is.character(model)) {
		index.mat<-matrix(0,2,k)
		
		if (model == "BM1"){
			np=1
			index<-matrix(TRUE,2,k)
			index.mat[1,1:k]<-np+1
			index.mat[2,1:k]<-1
			param.count<-np+1
			bool=TRUE
		}
		if (model == "BMS"){
			np=k
			index<-matrix(TRUE,2,k)
			index.mat[1,1:k]<-np+1
			index.mat[2,1:k]<-1:np
			param.count<-np+1
			bool=FALSE
		}
		if (model == "OU1"){
			np=2
			index<-matrix(TRUE,2,k)
			index.mat[1,1:k]<-1
			index.mat[2,1:k]<-2
			if(root.station==TRUE){
			param.count<-np+1
			}
			if(root.station==FALSE){
				param.count<-np+2
			}
			bool=root.station
		}
		if (model == "OUM"){
			np=2
			index<-matrix(TRUE,2,k)
			index.mat[1,1:k]<-1
			index.mat[2,1:k]<-2
			if(root.station==TRUE){
				param.count<-np+1
			}
			if(root.station==FALSE){
				param.count<-np+2
			}
			bool=root.station
		}
		if (model == "OUMV") {
			np=k+1
			index<-matrix(TRUE,2,k)
			index.mat[1,1:k]<-1
			index.mat[2,1:k]<-2:(k+1)
			if(root.station==TRUE){
				param.count<-np+1
			}
			if(root.station==FALSE){
				param.count<-np+2
			}
			bool=root.station
		}
		if (model == "OUMA") {
			np=k+1
			index<-matrix(TRUE,2,k)
			index.mat[1,1:k]<-1:k
			index.mat[2,1:k]<-k+1
			if(root.station==TRUE){
				param.count<-np+1
			}
			if(root.station==FALSE){
				param.count<-np+2
			}
			bool=root.station
		}
		if (model == "OUMVA") {
			np=k*2
			index<-matrix(TRUE,2,k)
			index.mat[index]<-1:(k*2)
			if(root.station==TRUE){
				param.count<-np+1
			}
			if(root.station==FALSE){
				param.count<-np+2
			}
			bool=root.station
		}
	}
	Rate.mat <- matrix(1, 2, k)
	#Likelihood function for estimating model parameters
	dev<-function(p){
		
		Rate.mat[] <- c(p, 1e-10)[index.mat]
		N<-length(x[,1])
		V<-varcov.ou(phy, edges, Rate.mat, root.state=root.state, simmap.tree=simmap.tree)
		W<-weight.mat(phy, edges, Rate.mat, root.state=root.state, simmap.tree=simmap.tree, assume.station=bool)
		
		theta<-pseudoinverse(t(W)%*%pseudoinverse(V)%*%W)%*%t(W)%*%pseudoinverse(V)%*%x
		
		DET<-determinant(V, logarithm=TRUE)

		logl<--.5*(t(W%*%theta-x)%*%pseudoinverse(V)%*%(W%*%theta-x))-.5*as.numeric(DET$modulus)-.5*(N*log(2*pi))
		
		return(-logl)
	}
	
	if(quiet==FALSE){
		cat("Begin subplex optimization routine -- Starting value:",ip, "\n")
	}
	
	lower = rep(lb, np)
	upper = rep(ub, np)
	
	opts <- list("algorithm"="NLOPT_LN_SBPLX", "maxeval"="1000000", "ftol_rel"=.Machine$double.eps^0.5)
	out = nloptr(x0=rep(ip, length.out = np), eval_f=dev, lb=lower, ub=upper, opts=opts)
	
	loglik <- -out$objective
	
	#Takes estimated parameters from dev and calculates theta for each regime:
	dev.theta<-function(p){
		
		tmp<-NULL
		Rate.mat[] <- c(p, 1e-10)[index.mat]
		
		N<-length(x[,1])
		V<-varcov.ou(phy, edges, Rate.mat, root.state=root.state, simmap.tree=simmap.tree)
		W<-weight.mat(phy, edges, Rate.mat, root.state=root.state, simmap.tree=simmap.tree, assume.station=bool)
		
		theta<-pseudoinverse(t(W)%*%pseudoinverse(V)%*%W)%*%t(W)%*%pseudoinverse(V)%*%x
		#Calculates the hat matrix:
		#H.mat<-W%*%pseudoinverse(t(W)%*%pseudoinverse(V)%*%W)%*%t(W)%*%pseudoinverse(V)
		#Standard error of theta -- uses pseudoinverse to overcome singularity issues
		se<-sqrt(diag(pseudoinverse(t(W)%*%pseudoinverse(V)%*%W)))
		#If plot.resid=TRUE, then the residuals are plotted with each regime designated -- colors start at 1 and go to n, where n is the total number of regimes
		tmp$res<-W%*%theta-x
		#Joins the vector of thetas with the vector of standard errors into a 2 column matrix for easy extraction at the summary stage
		tmp$theta.est<-cbind(theta,se)
		#Returns final GLS solution
		tmp
	}
	#Informs the user that the summarization has begun, output model-dependent and dependent on whether the root theta is to be estimated
	if(quiet==FALSE){
		cat("Finished. Summarizing results.", "\n")	
	}
	theta <- dev.theta(out$solution)
	#Calculates the Hessian for use in calculating standard errors and whether the maximum likelihood solution was found
	h <- hessian(x=out$solution, func=dev)
	#Using the corpcor package here to overcome possible NAs with calculating the SE
	solution<-matrix(out$solution[index.mat], dim(index.mat))
	solution.se<-matrix(sqrt(diag(pseudoinverse(h)))[index.mat], dim(index.mat))
	rownames(solution) <- rownames(solution.se) <- c("alpha","sigma.sq")
	if(simmap.tree==FALSE){
		colnames(solution) <- colnames(solution.se) <- levels(tot.states)
	}
	if(simmap.tree==TRUE){
		colnames(solution) <- colnames(solution.se) <- colnames(phy$mapped.edge)
	}				
	#Eigendecomposition of the Hessian to assess reliability of likelihood estimates
	hess.eig<-eigen(h,symmetric=TRUE)
	#If eigenvect is TRUE then the eigenvector and index matrix will appear in the list of objects 
	eigval<-signif(hess.eig$values,2)
	eigvect<-round(hess.eig$vectors, 2)
	obj = list(loglik = loglik, AIC = -2*loglik+2*param.count,AICc=-2*loglik+(2*param.count*(ntips/(ntips-param.count-1))),model=model,solution=solution, theta=theta$theta.est, solution.se=solution.se, tot.states=tot.states, index.mat=index.mat, simmap.tree=simmap.tree, opts=opts, data=data, phy=phy, root.station=root.station, lb=lower, ub=upper, iterations=out$iterations, res=theta$res, eigval=eigval, eigvect=eigvect) 
	class(obj)<-"OUwie"		
	return(obj)
}

print.OUwie<-function(x, ...){
	
	ntips=Ntip(x$phy)
	output<-data.frame(x$loglik,x$AIC,x$AICc,x$model,ntips, row.names="")
	names(output)<-c("-lnL","AIC","AICc","model","ntax")
	cat("\nFit\n")
	print(output)
	cat("\n")
					 
	if (is.character(x$model)) {
		if (x$model == "BM1" | x$model == "BMS"){
			param.est <- x$solution
			theta.mat <- matrix(t(x$theta[1,]), 2, length(levels(x$tot.states)))
			rownames(theta.mat)<-c("estimate", "se")
			if(x$simmap.tree==FALSE){
				colnames(theta.mat) <- levels(x$tot.states)
			}
			if(x$simmap.tree==TRUE){
				colnames(theta.mat) <- c(colnames(x$phy$mapped.edge))
			}
			cat("Rates\n")
			print(param.est)
			cat("\n")
			cat("Optima\n")
			print(theta.mat)
			cat("\n")
		}
		if (x$root.station == TRUE | x$root.station==FALSE){
			if (x$model == "OU1"){
				param.est<- x$solution
				theta.mat <- matrix(t(x$theta[1,]), 2, length(levels(x$tot.states)))
				rownames(theta.mat)<-c("estimate", "se")
				if(x$simmap.tree==FALSE){
					colnames(theta.mat) <- levels(x$tot.states)
				}
				if(x$simmap.tree==TRUE){
					colnames(theta.mat) <- c(colnames(x$phy$mapped.edge))
				}
				cat("Rates\n")
				print(param.est)
				cat("\n")
				cat("\nOptima\n")
				print(theta.mat)
				cat("\n")
			}
		}
		if (x$root.station == TRUE){
			if (x$model == "OUM"| x$model == "OUMV"| x$model == "OUMA" | x$model == "OUMVA"){
				param.est<- x$solution
				theta.mat<-matrix(t(x$theta), 2, length(levels(x$tot.states)))
				rownames(theta.mat)<-c("estimate", "se")
				if(x$simmap.tree==FALSE){
					colnames(theta.mat)<- levels(x$tot.states)
				}
				if(x$simmap.tree==TRUE){
					colnames(theta.mat) <- c(colnames(x$phy$mapped.edge))
				}
				cat("\nRates\n")
				print(param.est)
				cat("\n")
				cat("Optima\n")
				print(theta.mat)
				cat("\n")
			}
		}
		if (x$root.station == FALSE){
			if (x$model == "OUM"| x$model == "OUMV"| x$model == "OUMA" | x$model == "OUMVA"){ 
				print(x$theta)
				param.est<- x$solution
				theta.mat<-matrix(t(x$theta), 2, length(levels(x$tot.states))+1)
				rownames(theta.mat)<-c("estimate", "se")
				if(x$simmap.tree==FALSE){
					colnames(theta.mat)<-c("Root", levels(x$tot.states))
				}
				if(x$simmap.tree==TRUE){
					colnames(theta.mat)<-c("Root", colnames(phy$mapped.edge))
				}
				cat("\nRates\n")
				print(param.est)
				cat("\n")
				cat("Optima\n")
				print(theta.mat)
				cat("\n")
			}
		}		
	}
	if(any(x$eigval<0)){
		index.matrix <- x$index.mat
		rownames(index.matrix)<-c("alpha","sigma.sq")
		if(x$simmap.tree==FALSE){
			colnames(index.matrix) <- levels(x$tot.states)
		}
		if(x$simmap.tree==TRUE){
			colnames(index.matrix) <- colnames(x$phy$mapped.edge)
		}				
		#If any eigenvalue is less than 0 then the solution is not the maximum likelihood solution
		if (any(x$eigval<0)) {
			cat("The objective function may be at a saddle point -- check eigenvectors or try a simpler model", "\n")
		}
	}
	else{
		cat("Arrived at a reliable solution","\n")
	}
	
}
