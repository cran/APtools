#' @importFrom survival survfit
#' @importFrom survival coxph
#' @importFrom stats rexp
#' @importFrom stats quantile
#' @importFrom utils write.csv
#' @importFrom graphics par
#' @importFrom graphics plot
#' @importFrom graphics lines
#' @importFrom graphics legend
#' @export CompareAP

CompareAP <- function(status,marker1,marker2,stime=NULL,t0.list=NULL,method="perturbation",alpha=0.95,B=1000,weight=NULL,Plot=TRUE)
{
	############Checking the Formation############
	
	if(is.null(stime)&(!is.null(t0.list))){
		stop("When stime is NULL, t0.list should be NULL!\n") 
	}
	
	if(is.null(stime)){
		if((length(status)!=length(marker1))|(length(marker1)!=length(marker2))){
			stop("The length of each data is not equal!\n") 
		}
		data0=cbind(status,marker1,marker2)
		nn<-nrow(data0)
		vk = rep(1,nn)
		auc=ap=array(0,dim=c(B+1,2))
	
		if(method=="perturbation"){	
			vk1<-matrix(rexp(nn*B,1),nrow=nn,ncol=B)
			for(i in 1:2){
				dk=data0[,1];zk=data0[,i+1]			
				#auc[1,i] = sum((0.5*sum.I(zk,"<",zk,1*(dk==1)*vk)+0.5*sum.I(zk,"<=",zk,1*(dk==1)*vk))*(dk==0)*vk)/(sum(vk*(dk==1))*sum(vk*(dk==0))) 
				ap[1,i] =  sum(vk*(dk==1)*sum.I(zk,"<=",zk,1*(dk==1)*vk)/sum.I(zk,"<=",zk,vk),na.rm=T)/sum((dk==1)*vk) 
				####save true value
				cat("perturbation starts\n")
				#auc[2:(B+1),i] = apply((0.5*sum.I(zk,"<",zk,1*(dk==1)*vk1)+0.5*sum.I(zk,"<=",zk,1*(dk==1)*vk1))*(dk==0)*vk1,2,sum,na.rm=T)/(apply(vk1*(dk==1),2,sum)*apply(vk1*(dk==0),2,sum))
				ap[2:(B+1),i] = apply(vk1*(dk==1)*sum.I(zk,"<=",zk,1*(dk==1)*vk1)/sum.I(zk,"<=",zk,vk1),2,sum,na.rm=T)/apply((dk==1)*vk1,2,sum)	
			}
		}
	
		if(method=="bootstrap"){
			data_resam=array(0,dim=c(nn,ncol(data0),B+1))
			data_resam[,,1]=as.matrix(data0)
			for(k in 2:(B+1)){   
				index=sample(c(1:nn),nn,replace=TRUE)
				data_resam[,,k]=as.matrix(data0[as.vector(index),])
			} 
			for(i in 1:2){
				cat("Bootstrap starts\n")
				for(k in 1:(B+1)){
					dk=data_resam[,1,k];zk=data_resam[,i+1,k]
					#auc[k,i] = sum((0.5*sum.I(zk,"<",zk,1*(dk==1)*vk)+0.5*sum.I(zk,"<=",zk,1*(dk==1)*vk))*(dk==0)*vk)/(sum(vk*(dk==1))*sum(vk*(dk==0))) 
					ap[k,i] =  sum(vk*(dk==1)*sum.I(zk,"<=",zk,1*(dk==1)*vk)/sum.I(zk,"<=",zk,vk),na.rm=T)/sum((dk==1)*vk)  
				}
			}
		}
		event_rate=sum(status)/nn
		event_rate_sd=sqrt(event_rate*(1-event_rate)/nn)
		
		dap_summary=array(0,dim=c(5,3))
		
		######summary of propertion of cases
		dap_summary[1,1]<-event_rate
		dap_summary[1,2]<-event_rate-1.96*event_rate_sd
		dap_summary[1,3]<-event_rate+1.96*event_rate_sd
		######summary of AP
		for (i in 1:2){    
			dap_summary[i+1,1]<-ap[1,i]
			dap_summary[i+1,2]<-max(quantile(ap[,i],(1-alpha)/2,na.rm=T),0)
			dap_summary[i+1,3]<-min(quantile(ap[,i],(1+alpha)/2,na.rm=T),1)
		}
		#####summary of AP1-AP2#####
		if((mean(ap[,1])-mean(ap[,2]))>=0){
			dap1=ap[,1]-ap[,2]
			flag1="AP1-AP2"
		}
		else{
			dap1=ap[,2]-ap[,1]
			flag1="AP2-AP1"
		} 
		dap_summary[4,1]<-dap1[1]
		dap_summary[4,2]<-quantile(dap1,(1-alpha)/2,na.rm=T)
		dap_summary[4,3]<-quantile(dap1,(1+alpha)/2,na.rm=T)
		
		#####summary of AP1/AP2#####		
		if((mean(ap[,1])/mean(ap[,2]))>=1){
			dap2=ap[,1]/ap[,2]
			flag2="AP1/AP2"
		}
		else{
			dap2=ap[,2]/ap[,1]
			flag2="AP2/AP1"
		}
		dap_summary[5,1]<-dap2[1]
		dap_summary[5,2]<-quantile(dap2,(1-alpha)/2,na.rm=T)
		dap_summary[5,3]<-quantile(dap2,(1+alpha)/2,na.rm=T)
	
		colnames(dap_summary)<-c("point estimation",paste("Lower Limit(a=",alpha,")",sep=""),paste("Upper Limit(a=",alpha,")",sep=""))
		rownames(dap_summary)<-c("propertion of cases","AP1","AP2",flag1,flag2)
		write.csv(signif(dap_summary,3),file=paste("CompareAP_Binary_dap_summary(","method=",method,",B=",B,").csv",sep=""))
		return(list(dap_summary=signif(dap_summary,3)))
	}
	
	if(!is.null(stime)){
		if((length(stime)!=length(status))|(length(status)!=length(marker1))|(length(marker1)!=length(marker2))){
			stop("The length of each data is not equal!\n") 
		}
		if(is.null(t0.list)){
			stop("Please entry t0.list: prediction time intervals of interest for event time outcome!\n") 
		}				
		fit1=coxph(Surv(stime,status)~1)
		dfit1=survfit(fit1)
		tt=dfit1$time
		if(max(t0.list)>=max(tt)){
			stop("The prediction time intervals of interest are out of range!\n") 
		}	
		
		data0=cbind(stime,status,marker1,marker2)
		N_j=length(t0.list)
		nn<-nrow(data0)
		auc=ap=array(0,dim=c(B+1,N_j,2))
		Ti = data0[,1]; Di = data0[,2]
		vk = rep(1,nn)
		
		if(method=="perturbation"){
			for (j in 1:N_j){
				t0<-t0.list[j]		
				vk1<-matrix(rexp(nn*B,1),nrow=nn,ncol=B)
				if(is.null(weight)){		
				############Calculate the Weight############
					tt = c(t0,Ti[Ti<=t0])
					Wi = rep(0,length(Ti)); Vi=rep(1,length(Ti))
					tmpind = rank(tt)
					Ghat.tt = summary(survfit(Surv(Ti,1-Di)~1, se.fit=F, type='fl', weights=Vi), sort(tt))$surv[tmpind]
					Wi[Ti <= t0] = Di[Ti<=t0]/Ghat.tt[-1]; Wi[Ti >  t0] = 1/Ghat.tt[1]
					wk = Wi
					wk1=array(wk,dim=c(length(wk),B))
				}else{
					wk=weight
					wk1=array(wk,dim=c(length(wk),B))
				}
				for(i in 1:2){
					xk <- data0[,1]; zk <- data0[,i+2]; 
					#auc[1,j,i] = sum((0.5*sum.I(zk,"<=",zk,1*(xk<=t0)*wk*vk)+0.5*sum.I(zk,"<",zk,1*(xk<=t0)*wk*vk))*(xk>t0)*wk*vk)/(sum(vk*wk*(xk<=t0))*sum(vk*wk*(xk>t0))) 
					ap[1,j,i] = sum(wk*vk*(xk<=t0)*sum.I(zk,"<=",zk,1*(xk<=t0)*vk*wk)/sum.I(zk,"<=",zk,vk),na.rm=T)/sum((xk<=t0)*vk*wk)
					cat("perturbation starts\n")	
					#auc[2:(B+1),j,i]= apply(0.5*sum.I(zk,"<=",zk,1*(xk<=t0)*wk1*vk1)*(xk>t0)*wk1*vk1+0.5*sum.I(zk,"<",zk,1*(xk<=t0)*wk1*vk1)*(xk>t0)*wk1*vk1,2,sum,na.rm=T)/(apply(vk1*wk1*(xk<=t0),2,sum)*apply(vk1*wk1*(xk>t0),2,sum))
					ap[2:(B+1),j,i] = apply(wk1*vk1*(xk<=t0)*sum.I(zk,"<=",zk,1*(xk<=t0)*vk1*wk1)/sum.I(zk,"<=",zk,vk1),2,sum,na.rm=T)/apply((xk<=t0)*vk1*wk1,2,sum)
				}
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
			for (j in 1:N_j){
				t0<-t0.list[j]
				if(is.null(weight)){
					tt = c(t0,Ti[Ti<=t0])
					Wi = rep(0,length(Ti)); Vi=rep(1,length(Ti))
					tmpind = rank(tt)
					Ghat.tt = summary(survfit(Surv(Ti,1-Di)~1, se.fit=F, type='fl', weights=Vi), sort(tt))$surv[tmpind]
					Wi[Ti <= t0] = Di[Ti<=t0]/Ghat.tt[-1]; Wi[Ti >  t0] = 1/Ghat.tt[1]
					wkc = Wi
				}else{
					wkc=weight
				}
				for(i in 1:2){
					cat("Bootstrap starts\n")
					for(k in 1:(B+1)){
						wk=wkc[index[,k]]
						xk <- data_resam[,1,k]; zk <- data_resam[,i+2,k]; 
						#auc[k,j,i]= sum((0.5*sum.I(zk,"<=",zk,1*(xk<=t0)*wk*vk)+0.5*sum.I(zk,"<",zk,1*(xk<=t0)*wk*vk))*(xk>t0)*wk*vk)/(sum(vk*wk*(xk<=t0))*sum(vk*wk*(xk>t0))) 
						ap[k,j,i] = sum(wk*vk*(xk<=t0)*sum.I(zk,"<=",zk,1*(xk<=t0)*vk*wk)/sum.I(zk,"<=",zk,vk),na.rm=T)/sum((xk<=t0)*vk*wk)
					}
				}
			}
		}
		#use survival to find the corresponding event rate r based on t0
		fit1<-coxph(Surv(data0[,1],data0[,2])~1)
		dfit1<-survfit(fit1)
		tt<-dfit1$time
		rr <- 1-dfit1$surv
		pi.list <- t0.list*0
		for (l in 1:length(t0.list)){
			pi.list[l] = rr[which(tt>=t0.list[l])[1]]
		}

		dap_summary=matrix(0,nrow=N_j,ncol=14)
		dap_summary[,1]=t0.list
		dap_summary[,2]=pi.list
		######summary of AP
		for (j in 1:N_j){      
			dap_summary[j,3]<-ap[1,j,1]
			dap_summary[j,4]<-max(quantile(ap[,j,1],(1-alpha)/2,na.rm=T),0)
			dap_summary[j,5]<-min(quantile(ap[,j,1],(1+alpha)/2,na.rm=T),1)
			dap_summary[j,6]<-ap[1,j,2]
			dap_summary[j,7]<-max(quantile(ap[,j,2],(1-alpha)/2,na.rm=T),0)
			dap_summary[j,8]<-min(quantile(ap[,j,2],(1+alpha)/2,na.rm=T),1)
		}
		#####summary of AP1-AP2#####
		if((mean(ap[,,1])-mean(ap[,,2]))>=0){
			dap1=ap[,,1]-ap[,,2]
			flag1="AP1(t)-AP2(t)"
		}else{
			dap1=ap[,,2]-ap[,,1]
			flag1="AP2(t)-AP1(t)"
		}
		for (j in 1:N_j){ 
			dap_summary[j,9]<-dap1[1,j]
			dap_summary[j,10]<-quantile(dap1[,j],(1-alpha)/2,na.rm=T)
			dap_summary[j,11]<-quantile(dap1[,j],(1+alpha)/2,na.rm=T)
		}		
		#####summary of AP1/AP2#####		
		if((mean(ap[,,1])/mean(ap[,,2]))>=1){
			dap2=ap[,,1]/ap[,,2]
			flag2="AP1(t)/AP2(t)"
		}else{
			dap2=ap[,,2]/ap[,,1]
			flag2="AP2(t)/AP1(t)"
		}
		for (j in 1:N_j){
			dap_summary[j,12]<-dap2[1,j]
			dap_summary[j,13]<-quantile(dap2[,j],(1-alpha)/2,na.rm=T)
			dap_summary[j,14]<-quantile(dap2[,j],(1+alpha)/2,na.rm=T)
		}
		colnames(dap_summary)<-c("t0=","event rate","AP1(t)","(L,","U)","AP2(t)","(L,","U)",flag1,"(L,","U)",flag2,"(L,","U)")
		#rownames(dap_summary) = c(t0.list)
		write.csv(signif(dap_summary,3),file=paste("CompareAP_Survival_dap_summary(","method=",method,",B=",B,").csv",sep=""))

	
		###########plot##################
		if(Plot==TRUE){
			t0_l=seq(from=min(stime),to=max(stime),length.out=102)[c(-1,-102)]
			ap_plot=auc_plot=matrix(0,nrow=length(t0_l),ncol=2)
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
				xk=stime;zk=marker1
				ap_plot[j,1] = sum(wk*vk*(xk<=t0)*sum.I(zk,"<=",zk,1*(xk<=t0)*vk*wk)/sum.I(zk,"<=",zk,vk),na.rm=T)/sum((xk<=t0)*vk*wk)
				auc_plot[j,1] = sum((0.5*sum.I(zk,"<=",zk,1*(xk<=t0)*wk*vk)+0.5*sum.I(zk,"<",zk,1*(xk<=t0)*wk*vk))*(xk>t0)*wk*vk)/(sum(vk*wk*(xk<=t0))*sum(vk*wk*(xk>t0))) 
				zk=marker2
				ap_plot[j,2] = sum(wk*vk*(xk<=t0)*sum.I(zk,"<=",zk,1*(xk<=t0)*vk*wk)/sum.I(zk,"<=",zk,vk),na.rm=T)/sum((xk<=t0)*vk*wk)
				auc_plot[j,2] = sum((0.5*sum.I(zk,"<=",zk,1*(xk<=t0)*wk*vk)+0.5*sum.I(zk,"<",zk,1*(xk<=t0)*wk*vk))*(xk>t0)*wk*vk)/(sum(vk*wk*(xk<=t0))*sum(vk*wk*(xk>t0))) 
			}
			if((mean(ap_plot[,1])/mean(ap_plot[,2]))>=1){
				dap_plot=ap_plot[,1]/ap_plot[,2]
				flag2="AP1/AP2"
			}else{
				dap_plot=ap_plot[,2]/ap_plot[,1]
				flag2="AP2/AP1"
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
			par(mfrow=c(1,3))

			plot(t0_l,rep(0.5,length(t0_l)),type="l",xlim=c(0,max(t0_l)),ylim=c(0.5,max(auc_plot)),col="purple",lwd=2,xlab="Time",ylab="AUC",main="AUC vs t0",cex.main=1.5,cex.lab=1.2)
			lines(t0_l,auc_plot[,1],col="black",lwd=2)
			lines(t0_l,auc_plot[,2],col="red",lwd=2)
			legend("left",c("random marker","marker1","marker2"),bty="n",col=c("purple","black","red"),lwd=2,cex=1.2)  

			plot(t0_l,pi_l,type="l",xlim=c(0,max(t0_l)),ylim=c(0,max(ap_plot)),col="purple",lwd=2,xlab="Time",ylab="AP",main="AP vs t0",cex.main=1.5,cex.lab=1.2)
			lines(t0_l,ap_plot[,1],col="black",lwd=2)
			lines(t0_l,ap_plot[,2],col="red",lwd=2)
			legend("topleft",c("random marker","marker1","marker2"),bty="n",col=c("purple","black","red"),lwd=2,cex=1.2)  

			plot(t0_l,dap_plot,type="l",xlim=c(0,max(t0_l)),ylim=c(0,max(dap_plot)),lty=1,col="black",lwd=2,xlab="Time",ylab=flag2,main=paste(flag2,"vs t0"),cex.main=1.5,cex.lab=1.2)
						
		}

		return(list(dap_summary=signif(dap_summary,3)))
	}
  
}