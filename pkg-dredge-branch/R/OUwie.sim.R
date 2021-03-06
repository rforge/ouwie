##OUwie Simulator##

#written by Jeremy M. Beaulieu and Brian C. OMeara

#Simulates the Hansen model of continuous characters evolving under discrete selective 
#regimes. The input is a tree of class "phylo" that has the regimes as internal node labels 
#and a data file the contains the regime states for each species. The trait file must be in 
#the following order: Species names then Regime. The user must specify the parameters values
#for each simulation (i.e. alpha, sigma.sq, theta0, theta). 

##The following examples assume 2 selective regimes and different models can be specified:
#single rate Brownian motion BM1: alpha=c(0,0); sigma.sq=c(0.9); theta0=0; theta=c(0,0)
#two rate Brownian motion BMS: alpha=c(0,0); sigma.sq=c(0.45,.9); theta0=0; theta=c(0,0)
#global OU (OU1): alpha=c(0.1,0.1); sigma.sq=c(0.9,0.9); theta0=1; theta=c(1,1)
#normal OU (OUSM): alpha=c(0.1,0.1); sigma.sq=c(0.9,0.9); theta0=0; theta=c(1,2)
#multiple sigmas (OUSMV): alpha=c(0.1,0.1); sigma.sq=c(0.45,0.9); 
#multiple alphas (OUSMA): alpha=c(0.5,0.1); sigma.sq=c(0.9,0.9); theta0=0; theta=c(1,2)
#multiple alphas and sigmas (OUSMVA): alpha=c(0.5,0.1); sigma.sq=c(0.45,0.9); theta0=0; theta=c(1,2)

OUwie.sim <- function(phy, data=NULL, simmap.tree=FALSE, alpha, sigma.sq, theta0, theta){

	if(simmap.tree==FALSE){
		#This is annoying, but the second column has to be in there twice otherwise, error.
		data<-data.frame(data[,2], data[,2], row.names=data[,1])
		data<-data[phy$tip.label,]
		
		n=max(phy$edge[,1])
		ntips=length(phy$tip.label)
		
		int.states<-factor(phy$node.label)
		phy$node.label=as.numeric(int.states)
		tip.states<-factor(data[,1])
		data[,1]<-as.numeric(tip.states)
		k<-length(levels(int.states))
		
		regime=matrix(rep(0,(n-1)*k), n-1, k)
		
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
		
		#Resort the edge matrix so that it looks like the original matrix order
		edges=edges[sort.list(edges[,1]),]
		
		oldregime=root.state
		oldtime=1
		
		alpha=alpha
		sigma=sqrt(sigma.sq)
		theta=theta
		
		n.cov=matrix(rep(0,n*n), n, n)
		nodecode=matrix(c(ntips+1,1,1),1,3)
		
		x <- matrix(0, n, 1)
		TIPS <- 1:ntips
		ROOT <- ntips + 1L
		x[ROOT,] <- theta0
		
		for(i in 1:length(edges[,1])){
			anc = edges[i, 2]
			desc = edges[i, 3]
			
			newregime=which(edges[i,5:(k+4)]==1)
			current=edges[i,4]
			if(anc%in%nodecode[,1]){
				start=which(nodecode[,1]==anc)
				oldtime=nodecode[start,2]
				oldregime=nodecode[start,3]
			}
			else{
				newrow=c(anc,newtime,oldregime)
				nodecode=rbind(nodecode,newrow)
				oldtime=newtime
			}
			if(oldregime==newregime){
				newtime=current
				x[edges[i,3],]=x[edges[i,2],]*exp(-alpha[oldregime]*(oldtime-newtime))+(theta[oldregime])*(1-exp(-alpha[oldregime]*(oldtime-newtime)))+sigma[oldregime]*rnorm(1,0,1)*sqrt((1-exp(-2*alpha[oldregime]*(oldtime-newtime)))/(2*alpha[oldregime]))
			}
			else{
				newtime=current+((oldtime-current)/2)
				epoch1=x[edges[i,2],]*exp(-alpha[oldregime]*(oldtime-newtime))+(theta[oldregime])*(1-exp(-alpha[oldregime]*(oldtime-newtime)))+sigma[oldregime]*rnorm(1,0,1)*sqrt((1-exp(-2*alpha[oldregime]*(oldtime-newtime)))/(2*alpha[oldregime]))
				oldtime=newtime
				newtime=current
				x[edges[i,3],]=epoch1*exp(-alpha[newregime]*(oldtime-newtime))+(theta[newregime])*(1-exp(-alpha[newregime]*(oldtime-newtime)))+sigma[newregime]*rnorm(1,0,1)*sqrt((1-exp(-2*alpha[newregime]*(oldtime-newtime)))/(2*alpha[newregime]))
			}
			oldregime=newregime
		}
		x<-round(x,8)
		
		sim.dat<-matrix(,ntips,3)
		sim.dat<-data.frame(sim.dat)
		
		sim.dat[,1]<-phy$tip.label
		sim.dat[,2]<-data[,1]
		sim.dat[,3]<-x[TIPS]

		colnames(sim.dat)<-c("Genus_species","Reg","X")		
	}
	if(simmap.tree==TRUE){
		n=max(phy$edge[,1])
		ntips=length(phy$tip.label)
		
		k=length(colnames(phy$mapped.edge))

		regimeindex<-colnames(phy$mapped.edge)
		##Begins the construction of the edges matrix -- similar to the ouch format##
		#Makes a vector of absolute times in proportion of the total length of the tree
		branch.lengths=rep(0,(n-1))
		branch.lengths[(ntips+1):(n-1)]=branching.times(phy)[-1]/max(branching.times(phy))
		
		#Obtain root state and internal node labels
		root.state<-phy$map[[1]]
		
		#New tree matrix to be used for subsetting regimes
		edges=cbind(c(1:(n-1)),phy$edge,phy$edge.length)
		edges=edges[sort.list(edges[,3]),]
		
		edges[,4]=branch.lengths
		edges[,4]<-1-edges[,4]
		#Resort the edge matrix so that it looks like the original matrix order
		edges=edges[sort.list(edges[,1]),]
		
		oldregime=root.state
		oldtime=0
		
		alpha=alpha
		sigma=sqrt(sigma.sq)
		theta=theta
		
		n.cov=matrix(rep(0,n*n), n, n)
		nodecode=matrix(c(ntips+1,0,1),1,3)
		
		x <- matrix(0, n, 1)
		TIPS <- 1:ntips
		ROOT <- ntips + 1L
		x[ROOT,] <- theta0
		
		for(i in 1:length(edges[,1])){
			
			anc = edges[i, 2]
			desc = edges[i, 3]
			
			currentmap<-phy$maps[[i]]
			current=edges[i,4]
			
			if(anc%in%nodecode[,1]){
				start=which(nodecode[,1]==anc)
				oldtime=nodecode[start,2]
				oldregime=nodecode[start,3]
			}
			else{
				newrow=c(anc,newtime,oldregime)
				nodecode=rbind(nodecode,newrow)
				oldtime=newtime
			}
			
			for (regimeindex in 1:length(currentmap)){
				regimeduration<-currentmap[regimeindex]
				newtime<-oldtime+regimeduration
				regimenumber<-which(colnames(phy$mapped.edge)==names(currentmap)[regimeindex])
				x[edges[i,3],]=x[edges[i,3],]+x[edges[i,2],]*exp(-alpha[regimenumber]*(newtime-oldtime))+(theta[regimenumber])*(1-exp(-alpha[regimenumber]*(newtime-oldtime)))+sigma[regimenumber]*rnorm(1,0,1)*sqrt((1-exp(-2*alpha[regimenumber]*(newtime-oldtime)))/(2*alpha[regimenumber]))
				oldtime<-newtime
				newregime<-regimenumber
			}
			oldregime=newregime
		}
		x<-round(x,8)
		
		sim.dat<-matrix(,ntips,2)
		sim.dat<-data.frame(sim.dat)
		
		sim.dat[,1]<-phy$tip.label
		sim.dat[,2]<-x[TIPS]
		
		colnames(sim.dat)<-c("Genus_species","X")
	}
	sim.dat
}



