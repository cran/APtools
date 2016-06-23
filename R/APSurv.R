#' @importFrom survival survfit
#' @importFrom survival coxph
#' @importFrom stats rexp
#' @importFrom stats quantile
#' @importFrom utils write.csv
#' @importFrom graphics plot
#' @importFrom graphics lines
#' @importFrom graphics legend
#' @export APSurv

APSurv <- function(stime,status,marker,t0.list,cut.values=NULL,method="perturbation",alpha=0.95,B=1000,weight=NULL,Plot=TRUE)
{	

	############Checking the Formation############
	
	if((length(stime)!=length(status))|(length(status)!=length(marker))){
		stop("The length of each data is not equal!\n") 
	}
	fit1=coxph(Surv(stime,status)~1)
	dfit1=survfit(fit1)
	tt=dfit1$time
	if(max(t0.list)>=max(tt)){
		stop("The prediction time intervals of interest are out of range!\n") 
	}	
	
	
	N_j=length(t0.list)
	data0=cbind(stime,status,marker)
	nn=nrow(data0)
	auc=ap=ap_event=array(0,dim=c(B+1,N_j))
	
	############Set Cut-Off Value############
	xk=Ti=stime;zk=marker;Di=status
	vk = rep(1,nn)
	if(!is.null(cut.values)){
		scl=cut.values
		PPV=TPF=array(0,dim=c(length(scl),N_j+1))
		PPV[,1]=TPF[,1]=scl
	}
	

	if(method=="perturbation"){
		for (j in 1:N_j){
			t0=t0.list[j]
			vk1=matrix(rexp(nn*B,1),nrow=nn,ncol=B)
			
			if(is.null(weight)){
			############Calculate the Weight############
				tt = c(t0,Ti[Ti<=t0])
				Wi = rep(0,length(Ti)); Vi=rep(1,length(Ti))
				tmpind = rank(tt)
				Ghat.tt = summary(survfit(Surv(Ti,1-Di)~1, se.fit=F, type='fl', weights=Vi), sort(tt))$surv[tmpind]
				Wi[Ti <= t0] = Di[Ti<=t0]/Ghat.tt[-1]; Wi[Ti >  t0] = 1/Ghat.tt[1]
				wk = Wi
				wk1 = array(wk,dim=c(length(wk),B))
			}else{
				wk = weight
				wk1 = array(wk,dim=c(length(wk),B))
			}
			############Point estimation############
			
			if(!is.null(cut.values)){
				St0.Fcl = sum.I(scl,">=",zk,1*(xk>t0)*wk*vk)/sum(wk*vk) 
				Ft0.Fcl = sum.I(scl,">=",zk,1*(xk<=t0)*wk*vk)/sum(wk*vk) ## hh
				Fcl = sum.I(scl,">=",zk)/nn # ss	
				St0 = max(St0.Fcl); Ft0 = 1-St0 ## St0 = P(T> t0); Ft0 = P(T<=t0) = pi
				PPV[,j+1] = (Ft0-Ft0.Fcl)/(1-Fcl); if (sum(Fcl==1)>0) PPV.cl[Fcl==1] = NA  ## P(T<=t0|Y> cl)
				TPF[,j+1] = (Ft0-Ft0.Fcl)/Ft0     ## P(Y> cl|T<=t0)  
			}
		
			auc[1,j] = sum((0.5*sum.I(zk,"<=",zk,1*(xk<=t0)*wk*vk)+0.5*sum.I(zk,"<",zk,1*(xk<=t0)*wk*vk))*(xk>t0)*wk*vk)/(sum(vk*wk*(xk<=t0))*sum(vk*wk*(xk>t0))) 
			ap[1,j] = sum(wk*vk*(xk<=t0)*sum.I(zk,"<=",zk,1*(xk<=t0)*vk*wk)/sum.I(zk,"<=",zk,vk),na.rm=T)/sum((xk<=t0)*vk*wk)
			ap_event[1,j] = ap[1,j]*nn/sum((xk<=t0)*vk*wk)
			# ap_event[1,j] = ap[1,j]/pi.list[j]
			
			
			############Perturbation Process############
			cat("perturbation starts\n")
			
			auc[2:(B+1),j]= apply(0.5*sum.I(zk,"<=",zk,1*(xk<=t0)*wk1*vk1)*(xk>t0)*wk1*vk1+0.5*sum.I(zk,"<",zk,1*(xk<=t0)*wk1*vk1)*(xk>t0)*wk1*vk1,2,sum,na.rm=T)/(apply(vk1*wk1*(xk<=t0),2,sum)*apply(vk1*wk1*(xk>t0),2,sum))
			ap[2:(B+1),j] = apply(wk1*vk1*(xk<=t0)*sum.I(zk,"<=",zk,1*(xk<=t0)*vk1*wk1)/sum.I(zk,"<=",zk,vk1),2,sum,na.rm=T)/apply((xk<=t0)*vk1*wk1,2,sum)
			ap_event[2:(B+1),j] = ap[2:(B+1),j]*nn/sum((xk<=t0)*vk*wk)
			# ap_event[2:(B+1),j] = ap[2:(B+1),j]/pi.list[j]
		}	
	}
	
	if(method=="bootstrap"){
		data_resam=array(0,dim=c(nn,ncol(data0),B+1))
		data_resam[,,1]=as.matrix(data0)
		index=matrix(0,nrow=nn,ncol=B+1)
		index[,1]=seq(from=1,to=nn,length=nn)
		for(k in 2:(B+1)){   
			index[,k]=sample(c(1:nn),nn,replace=TRUE)
			data_resam[,,k]=as.matrix(data0[as.vector(index[,k]),])		
		}
		pi.list_bs=matrix(0,nrow=B+1,ncol=N_j)
		
		for (j in 1:N_j){
			t0<-t0.list[j]
			if(is.null(weight)){
				tt = c(t0,Ti[Ti<=t0])
				Vi = rep(1,length(Ti));Wi = rep(0,length(Ti));tmpind = rank(tt)
				Ghat.tt = summary(survfit(Surv(Ti,1-Di)~1, se.fit=F, type='fl', weights=Vi), sort(tt))$surv[tmpind]
				Wi[Ti <= t0] = Di[Ti<=t0]/Ghat.tt[-1]; Wi[Ti >  t0] = 1/Ghat.tt[1]
				wkc = Wi
			}else{
				wkc=weight
			}
			if(!is.null(cut.values)){
				wk=wkc
				St0.Fcl = sum.I(scl,">=",zk,1*(xk>t0)*wk*vk)/sum(wk*vk) 
				Ft0.Fcl = sum.I(scl,">=",zk,1*(xk<=t0)*wk*vk)/sum(wk*vk) ## hh
				Fcl = sum.I(scl,">=",zk)/nn # ss	
				St0 = max(St0.Fcl); Ft0 = 1-St0 ## St0 = P(T> t0); Ft0 = P(T<=t0) = pi
				PPV[,j+1] = (Ft0-Ft0.Fcl)/(1-Fcl); if (sum(Fcl==1)>0) PPV.cl[Fcl==1] = NA  ## P(T<=t0|Y> cl)
				TPF[,j+1] = (Ft0-Ft0.Fcl)/Ft0     ## P(Y> cl|T<=t0)  
			}
			
			cat("Bootstrap starts\n")
			for(k in 1:(B+1)){
				wk=wkc[index[,k]]
				xk <- data_resam[,1,k]; zk <- data_resam[,3,k]; dk <- data_resam[,2,k]
				auc[k,j]= sum((0.5*sum.I(zk,"<=",zk,1*(xk<=t0)*wk*vk)+0.5*sum.I(zk,"<",zk,1*(xk<=t0)*wk*vk))*(xk>t0)*wk*vk)/(sum(vk*wk*(xk<=t0))*sum(vk*wk*(xk>t0))) 
				ap[k,j] = sum(wk*vk*(xk<=t0)*sum.I(zk,"<=",zk,1*(xk<=t0)*vk*wk)/sum.I(zk,"<=",zk,vk),na.rm=T)/sum((xk<=t0)*vk*wk)
				ap_event[k,j] = ap[k,j]*nn/sum((xk<=t0)*vk*wk)
				# surv.out = survfit(Surv(xk,dk)~1)
				# tt = surv.out$time
				# rr = 1-surv.out$surv
				# pi3.out[k,j] = rr[which(tt>t0)[1]-1]
				# ap_event[k,j] = ap[k,j]/pi.list[k,j]					
			}
		}
	}

	####use survival to find the corresponding event rate r based on t0
	fit1<-coxph(Surv(data0[,1],data0[,2])~1)
	dfit1<-survfit(fit1)
	tt<-dfit1$time
	rr = 1-dfit1$surv
	pi.list = t0.list*0
	for (l in 1:length(t0.list)){
		pi.list[l] = rr[which(tt>=t0.list[l])[1]]
	}	
	
	auc_summary=array(0,dim=c(N_j,5))
	ap_summary=array(0,dim=c(N_j,8))
	auc_summary[,1]=ap_summary[,1]=t0.list
	auc_summary[,2]=ap_summary[,2]=signif(pi.list,3)
   
	#####summary of AUC#####
	for (j in 1:N_j){   
		auc_summary[j,3]=signif(auc[1,j],3)
		auc_summary[j,4]=signif(max(quantile(auc[,j],(1-alpha)/2,na.rm=T),0),3)
		auc_summary[j,5]=signif(min(quantile(auc[,j],(1+alpha)/2,na.rm=T),1),3)
	}

	######summary of AP
	for (j in 1:N_j){
		ap_summary[j,3]=signif(ap[1,j],3)
		ap_summary[j,4]=signif(max(quantile(ap[,j],(1-alpha)/2,na.rm=T),0),3)
		ap_summary[j,5]=signif(min(quantile(ap[,j],(1+alpha)/2,na.rm=T),1),3)
		ap_summary[j,6]=signif(ap_event[1,j],3)
		ap_summary[j,7]=signif(quantile(ap_event[,j],(1-alpha)/2,na.rm=T),3)
		ap_summary[j,8]=signif(quantile(ap_event[,j],(1+alpha)/2,na.rm=T),3)
	}
   
	colnames(auc_summary)=c("t0=","event rate","AUC(t)",paste("Lower Limit(a=",alpha,")",sep=""),paste("Upper Limit(a=",alpha,")",sep=""))
	#rownames(auc_summary) = c(t0.list)
	write.csv(auc_summary,file=paste("APSurv_auc_summary(","method=",method,",B=",B,").csv",sep=""))

	colnames(ap_summary)=c("t0=","event rate","AP(t)",paste("Lower Limit(a=",alpha,")",sep=""),paste("Upper Limit(a=",alpha,")",sep=""),"AP/(event rate)",paste("Lower Limit(a=",alpha,")",sep=""),paste("Upper Limit(a=",alpha,")",sep=""))
	#rownames(ap_summary) = c(t0.list)
	write.csv(ap_summary,file=paste("APSurv_ap_summary(","method=",method,",B=",B,").csv",sep=""))		
		
	###########plot##################
	if(Plot==TRUE){

		t0_l=seq(from=min(stime),to=max(stime),length.out=102)[c(-1,-102)]
		ap_plot=rep(0,length(t0_l))
		xk=stime;zk=marker
		for (j in 1:length(t0_l)){
			t0<-t0_l[j]
			if(is.null(weight)){
			############Calculate the Weight############
				tt = c(t0,Ti[Ti<=t0])
				Wi = rep(0,length(Ti)); Vi=rep(1,length(Ti))
				tmpind = rank(tt)
				Ghat.tt = summary(survfit(Surv(Ti,1-Di)~1, se.fit=F, type='fl', weights=Vi), sort(tt))$surv[tmpind]
				Wi[Ti <= t0] = Di[Ti<=t0]/Ghat.tt[-1]; Wi[Ti >  t0] = 1/Ghat.tt[1]
				wk = Wi
			}else{
				wk = weight
			}
			ap_plot[j]= sum(wk*vk*(xk<=t0)*sum.I(zk,"<=",zk,1*(xk<=t0)*vk*wk)/sum.I(zk,"<=",zk,vk),na.rm=T)/sum((xk<=t0)*vk*wk)
		}
		
		#use survival to find the corresponding event rate r based on t0
		fit1<-coxph(Surv(data0[,1],data0[,2])~1)
		dfit1<-survfit(fit1)
		tt<-dfit1$time
		rr <- 1-dfit1$surv
		pi_l <- t0_l*0
		for (l in 1:length(t0_l)){
			pi_l[l] = rr[which(tt>=t0_l[l])[1]]
		}
		
		###########plot##################
		plot(t0_l,pi_l,type="l",xlim=c(0,max(t0_l)),ylim=c(0,max(ap_plot)),col="purple",lwd=2,xlab="Time",ylab="AP",main="AP vs t0",)
		lines(t0_l,ap_plot,col="red",lwd=2)
		legend("topleft",c("random marker",colnames(data0)[3]),col=c("purple","red"),lwd=2,cex=1.2)  		
	}
	
	if(!is.null(cut.values)){
		colnames(PPV)=c("cut.off values",paste("t0=",c(t0.list)))
		colnames(TPF)=c("cut.off values",paste("t0=",c(t0.list)))
		write.csv(signif(PPV,3),file=paste("APSurv_PPV.csv",sep=""))
		write.csv(signif(TPF,3),file=paste("APSurv_TPF.csv",sep=""))
		return(list(PPV=signif(PPV,3),TPF=signif(TPF,3),ap_summary=ap_summary,auc_summary=auc_summary))
	}
	return(list(ap_summary=ap_summary,auc_summary=auc_summary))
}