#do mcmc on the discrete parameters of the model; use optimization on the other parameters. Saves us from having to do the rjmcmc, since all the parameters are always there.
library(coda)


#allow users to pass priors as vectors of the probability of 0, 1, 2, etc. up to the maximum number of, say, k.alpha
OUwie.dredge.mcmc<-function(phy,data, prior.k.theta, prior.k.sigma, prior.k.alpha, prior.nonmatching.theta, prior.nonmatching.sigma, prior.nonmatching.alpha, ngen=10000000, sample.freq=100, print.freq=10, summary.freq=10000, root.station=TRUE, ip=1, lb=0.000001, ub=1000,maxeval=500, samplesfile="mcmc.samples.txt",likelihoodfree=FALSE) {
	data2<-data.frame(as.character(data[,1]),sample(c(1:2),length(data[,1]), replace=T),data[,2],stringsAsFactors=FALSE)
	phy$node.label<-sample(c(1:2),phy$Nnode, replace=T)
	cat("Initializing...","\n")
	start<-OUwie(phy,data2,model=c("OU1"),plot.resid=FALSE, quiet=TRUE)
	nodes<-Nnode(phy,internal.only=FALSE)
	ip<-matrix(c(rep(start$theta[,1],nodes),rep(start$Param.est[2],nodes),rep(start$Param.est[1],nodes)),nrow=1,ncol=3*nodes) #BM1

	data<-data.frame(data[,2], data[,2], row.names=data[,1])
	data<-data[phy$tip.label,]
	bayes.individual<-c(rep(0,nodes),rep(0,nodes),rep(0,nodes)) #BM1
	bayes.individual[c(nodes,2*nodes,3*nodes)]<-c(1,1,-1)
	current.state<-measure.proposal(phy=phy, data=data, new.individual=bayes.individual, prior.k.theta=prior.k.theta, prior.k.sigma=prior.k.sigma, prior.k.alpha=prior.k.alpha, prior.nonmatching.theta=prior.nonmatching.theta, prior.nonmatching.sigma=prior.nonmatching.sigma, prior.nonmatching.alpha=prior.nonmatching.alpha, root.station=root.station, lb=lb, ub=ub, ip=ip, maxeval=maxeval,likelihoodfree=likelihoodfree,proposal.ratio=1)
	store.state(current.state,samplesfile,generation=0)
	cat("\nNote that this is loglik, not negative loglik. Higher is better (true for posterior, too)\n")
	labels<-c("generation","posterior","loglik", "prior", "k.theta", "k.sigma", "k.alpha", "nRegimes")
	cat("\n",labels)
	write.table(matrix(labels,nrow=1),file=samplesfile,quote=F,sep="\t",row.name=F,col.name=F,append=F)
	print.state(current.state,samplesfile,generation=0)

	for(generation in sequence(ngen)) {
		new.transform=transform.single(current.state$new.individual,phy)		
		next.state<-measure.proposal(phy=phy, data=data, new.individual=new.transform$new.individual, prior.k.theta=prior.k.theta, prior.k.sigma=prior.k.sigma, prior.k.alpha=prior.k.alpha, prior.nonmatching.theta=prior.nonmatching.theta, prior.nonmatching.sigma=prior.nonmatching.sigma, prior.nonmatching.alpha=prior.nonmatching.alpha, root.station=root.station, lb=lb, ub=ub, ip=ip, maxeval=maxeval,likelihoodfree=likelihoodfree,proposal.ratio=new.transform$proposal.ratio)
		#cat("\nproposed posterior and like and prior: ",next.state$posterior,next.state$loglik,next.state$prior)
		#cat("\nprevious posterior and like and prior: ",current.state$posterior,current.state$loglik,current.state$prior)
		#cat("\nposterior diff exp(new - current): ",exp(next.state$posterior - current.state$posterior))

		if (is.finite(next.state$posterior) && (exp(next.state$posterior - current.state$posterior) >= (runif(1,0,1)))) {
			current.state<-next.state
		}
		if (generation %% sample.freq == 0) {
			store.state(current.state,samplesfile,generation)
		}
		if (generation %% print.freq == 0) {
			print.state(current.state,samplesfile,generation)
		}
		if (generation %% summary.freq == 0) {
			summarize.results(samplesfile=samplesfile)
		}
	}
}

truncated.geometric.prior<-function(min=1,max=5,prob=0.5) {
	if(min<0) {
		stop("min in truncated.geometric.prior should be non-negative")
	}
	params<-c(min:max)
	probs<-sapply(params,dgeom,prob)
	probs<-probs/sum(probs) #remember this could be from min=2 to max=4. We want to start at zero
	allprobs<-rep(0,max+1) #first element is zero
	allprobs[params+1]<-probs #deals with the zero offset
	return(allprobs)
}

truncated.uniform.prior<-function(min=1,max=5) {
	if(min<0) {
		stop("min in truncated.uniform.prior should be non-negative")
	}
	params<-c(min:max)
	probs<-rep(1,max-min+1)
	probs<-probs/sum(probs)
	allprobs<-rep(0,max+1) #first element is zero
	allprobs[params+1]<-probs #deals with the zero offset
	return(allprobs)	
}

#transforms just a single cell
transform.single<-function(bayes.individual,phy) {
	valid.transform<-FALSE
	while(!valid.transform) { #cool thing about doing this is that the hastings ratio is one: prob of going from valid to invalid is zero, as is the reverse rate
		new.individual<-bayes.individual
		focal.param<-sample.int(length(new.individual),1)
		new.individual[focal.param]<-new.individual[focal.param]+sign(rnorm(1))
		valid.transform<-valid.individual(new.individual,phy)
	}
	return(list(new.individual=new.individual,proposal.ratio=1))
}



transform.clade<-function(bayes.individual,phy) {
	valid.transform<-FALSE
	while(!valid.transform) { #cool thing about doing this is that the hastings ratio is one: prob of going from valid to invalid is zero, as is the reverse rate
		new.individual<-bayes.individual
		
		#do things here
		
		
		valid.transform<-valid.individual(new.individual,phy)
	}
	return(list(new.individual=new.individual,proposal.ratio=1))
}

count.nonmatching<-function(bayes.individual) { #excludes root
	stepLength<-length(bayes.individual)/3
	bayes.individual<-bayes.individual[c(-stepLength, -2*stepLength, -3*stepLength)] #cutting out the root, as it can't be matching. So a nonzero prior of 0 nonmatching can match real data
	stepLength<-length(bayes.individual)/3
	results<-c(length(which(bayes.individual[1:stepLength] != 0)),  length(which(bayes.individual[(1+stepLength):(2*stepLength)] != 0)) , length(which(bayes.individual[((2*stepLength) + 1):(3*stepLength)] != 0)) )
	names(results)<-c("n.nonmatching.theta","n.nonmatching.sigma","n.nonmatching.alpha")
	return(results)
}

prior.prob<-function(bayes.individual,prior.k.theta,prior.k.sigma, prior.k.alpha, prior.nonmatching.theta, prior.nonmatching.sigma, prior.nonmatching.alpha, phy) {
	k.list<-param.count(bayes.individual, phy)
	k.number.prior<-vector.index.safe.offset(k.list$k.theta,prior.k.theta) * vector.index.safe.offset(k.list$k.sigma,prior.k.sigma) * vector.index.safe.offset(k.list$k.alpha,prior.k.alpha) 
	nonmatching<-count.nonmatching(bayes.individual)
	k.nonmatching.prior<-vector.index.safe.offset(nonmatching[1], prior.nonmatching.theta) * vector.index.safe.offset(nonmatching[2], prior.nonmatching.theta) * vector.index.safe.offset(nonmatching[3], prior.nonmatching.theta)
	return(k.number.prior * k.nonmatching.prior)
}

vector.index.safe<-function(x,vec) {
	if (x>length(vec)) {
		return(0)
	}
	else {
		return(vec[x])
	}
}

vector.index.safe.offset<-function(x,vec) {
	x<-x+1 #deals with the indexing by zero but R indexing by 1
	if (x>length(vec)) {
		return(0)
	}
	else {
		return(vec[x])
	}
}


measure.proposal<-function(phy, data, new.individual, prior.k.theta,prior.k.sigma, prior.k.alpha, prior.nonmatching.theta, prior.nonmatching.sigma, prior.nonmatching.alpha, maxeval, root.station, lb, ub, ip,likelihoodfree=FALSE,proposal.ratio) {
	edge.mat.all<-edge.mat(phy,new.individual)
	loglik=0.5
	if(!likelihoodfree) {
		loglik<-dev.optimize(edges.ouwie=edge.mat.all$edges.ouwie, regime.mat=edge.mat.all$regime, data=data,maxeval=maxeval, root.station=root.station,lb=lb, ub=ub, ip=ip, phy=phy)$loglik 
	}
	prior<-prior.prob(new.individual,prior.k.theta,prior.k.sigma, prior.k.alpha, prior.nonmatching.theta, prior.nonmatching.sigma, prior.nonmatching.alpha, phy)
	posterior<-loglik+log(prior)+log(proposal.ratio)
	k.list<-param.count(new.individual,phy)
	return(list(new.individual=new.individual, loglik=loglik, prior=prior, posterior=posterior, k.theta=k.list$k.theta, k.sigma=k.list$k.sigma, k.alpha=k.list$k.alpha,nRegimes=dim(edge.mat.all$regime)[1],proposal.ratio=proposal.ratio))
}

store.state<-function(current.state,samplesfile,generation) {
	tmp<-c(generation,current.state$posterior,current.state$loglik, current.state$prior, current.state$k.theta, current.state$k.sigma, current.state$k.alpha, current.state$nRegimes, current.state$new.individual)
	write.table(matrix(tmp,nrow=1),file=samplesfile,quote=F,sep="\t",row.name=F,col.name=F,append=T)
}

print.state<-function(current.state,samplesfile,generation) {
	cat("\n",c(generation,current.state$posterior,current.state$loglik, current.state$prior, current.state$k.theta, current.state$k.sigma, current.state$k.alpha, current.state$nRegimes))
}

summarize.results<-function(samplesfile="mcmc.samples.txt",q=0.4,r=0.4,verbose=TRUE) {
	cat("\n---------------------------------------------------------------------------------\nDoing summaries\n\nLoading samples...\n")
	samples<-read.table(file=samplesfile, skip=1) #skip the first line, as it has a header but not for every column
	names(samples)[1:8 ] <- c("generation","posterior","loglik","prior","k.theta","k.sigma","k.alpha","nRegimes")
	cat("\nCalculating effective sizes...\n")
	effectiveSizes<-effectiveSize(as.mcmc(samples))
	namedEffective<-effectiveSizes[2:8]
	namedEffective<-c(namedEffective,c(min(effectiveSizes[9:length(effectiveSizes)]), max(effectiveSizes[9:length(effectiveSizes)]), median(effectiveSizes[9:length(effectiveSizes)])))
	names(namedEffective)<-c(names(namedEffective)[1:7],"min ESS","max ESS","median ESS")
	cat("\nEstimating burnin...\n")
	rafteryResults<-raftery.diag(as.mcmc(samples[,5:dim(samples)[2]]),q=q,r=r)
	requiredGen<-max(rafteryResults$resmatrix[,2],na.rm=TRUE)
	burnin<-max(rafteryResults$resmatrix[,1],na.rm=TRUE)
	cat("\nCalculating estimated values...\n")
	estimates.noburn<-apply(samples[,2:8],MARGIN=2,quantile,probs=c(0,0.025,0.05,0.5,0.95,0.975,1),na.rm=TRUE)
	estimates.burn<-estimates.noburn
	estimates.burn<-NA
	try(estimates.burn<-apply(samples[burnin:(dim(samples)[1]),2:8],MARGIN=2,quantile,probs=c(0,0.025,0.05,0.5,0.95,0.975,1),na.rm=TRUE))
	if(verbose) {
		cat("\nEstimates no burnin\n")
		print(estimates.noburn)
		cat("\nEstimates with burnin\n")
		print(estimates.burn)
		cat("\nEffective sample sizes\n")
		print(namedEffective)
		cat("\n\nEstimate for required number of samples: ",requiredGen)
		cat("\n\nEstimate for required burnin: ",burnin, "(note that ",dim(samples)[1]," samples, representing ",max(samples[,1],na.rm=TRUE),"generations, have completed so far)")
	}
	cat("\n---------------------------------------------------------------------------------\n")
	return(list(estimates.noburn=estimates.noburn,estimates.burn=estimates.burn,effectiveSizes=effectiveSizes,namedEffective=namedEffective,requiredGen=requiredGen,burnin=burnin))
}