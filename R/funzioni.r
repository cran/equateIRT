
linkp<-function(coef)
{
	l<-length(coef)
	out<-matrix(NA,l,l)
	for(i in 1:l)
		for (j in 1:l)
			out[i,j]<-sum(rownames(coef[[i]])%in%rownames(coef[[j]]))
	return(out)
}


modIRT<-function(coef,var=NULL,names=NULL,ltparam=TRUE,lparam=TRUE,display=TRUE,digits=2)
{
	coef<-lapply(coef,FUN=function(x) x<-as.matrix(x[,colSums(x==1)!=nrow(x)]))
	coef<-lapply(coef,FUN=function(x) x<-as.matrix(x[,colSums(x==0)!=nrow(x)]))
	itmp<-sapply(coef,ncol)
	if (length(unique(itmp))!=1) stop("Mixed model types not allowed")
	itmp<-itmp[1]
	if (itmp==1 | itmp==2) lparam<-FALSE
	lc<-length(coef)
	if (is.null(names)) names<-paste("T",1:lc,sep="")
	mods<-list()
	for (i in 1:lc) {
		c0<-c1<-c2<-NULL
		D1<-D2<-matrix(NA,0,0)
		coefi<-coef[[i]]
		if (!is.null(var)) vari<-var[[i]]
		else vari<-NULL
		if (!is.null(vari)) if (any(is.na(vari))) vari<-NULL
		n1<-nrow(coefi)
		if (itmp==1)
			c1<-coefi[,1]
		if (itmp==2) { 
			c1<-coefi[,1]
			c2<-coefi[,2]
		}
		if (itmp==3) { 
			c0<-coefi[,1]
			c1<-coefi[,2]
			c2<-coefi[,3]
		}
		if (ltparam & itmp==1) {
			c1<--c1
		}
		if (ltparam & itmp>1) {
			D1<-rbind(diag(-1/c2),diag(c1/c2^2))
			D1<-cbind(D1,rbind(matrix(0,n1,n1),diag(1,n1)))
			c1<--c1/c2
		}
		if (lparam) {
			D2<-diag(exp(c0)/(1+exp(c0))^2)
			c0<-exp(c0)/(1+exp(c0))
		}
		if (!lparam & itmp==3) {
			D2<-diag(1,n1)
		}
		if (!is.null(vari)) {
			if ((ltparam | lparam) & itmp>1) {
				D1<-blockdiag(D2,D1)
				vari<-t(D1)%*%vari%*%D1
			}
		}
		if (itmp==3) names(c0)<-paste("Gussng",names(c0),sep=".")
		names(c1)<-paste("Dffclt",names(c1),sep=".")
		if (itmp>1) names(c2)<-paste("Dscrmn",names(c2),sep=".")
		coefi<-c(c0,c1,c2)
		if (!is.null(vari)) rownames(vari)<-colnames(vari)<-names(coefi)
		mods[[i]]<-list(coef=coefi,var=vari,itmp=itmp)
		if (display) {
			if (!is.null(var)) out<-cbind(coefi,diag(vari)^0.5)
			else out<-cbind(coefi,NA)
			out<-round(out,digits)
			colnames(out)<-c("value","std.err")
			cat("Form:",names[i],"\n")
			print(out)
		}
	}
	names(mods)<-names
	class(mods)<-"modIRT"
	return(mods)
}


irtp1<-function(ab,diff,discr,guess,D)
{
	elp<-exp(D*discr*(ab-diff))
	guess+(1-guess)*elp/(1+elp)
}



obj<-function(eqc,P2,ab,a1,b1,c1,met,itmp,wt,D=D){
	ifelse(itmp==1,A<-1,A<-eqc[1])
	B<-eqc[2]

	b12<-A*b1+B
	a12<-a1/A

	ni<-ncol(P2)
	P1<-matrix(NA,length(ab),ni)
	for (i in 1:ni)
		P1[,i]<-irtp1(ab,diff=b12[i],discr=a12[i],guess=c1[i],D=D)

	if (met=="Haebara") f<-0.5*sum(rowSums((P2-P1)^2)*wt)
	if (met=="Stocking-Lord") f<-0.5*sum(((rowSums(P2)-rowSums(P1))^2)*wt)
	return(f)
}


direc<-function(mod1,mod2,method="mean-mean",suff1=".1",suff2=".2",D=1,quadrature=TRUE,nq=30)
{
	if (method!="mean-mean" & method!="mean-sigma" & method!="mean-gmean" & method!="Haebara" & method!="Stocking-Lord") warning("Method not implemented.")
	name1<-names(mod1)
	name2<-names(mod2)
	forms<-paste(name1,name2,sep=".")
	mod1<-mod1[[1]]
	mod2<-mod2[[1]]
	tab1<-data.frame(value1=mod1$coef)
	tab2<-data.frame(value2=mod2$coef)
	var1<-mod1$var
	var2<-mod2$var
	itmp1<-mod1$itmp
	itmp2<-mod2$itmp
	if (itmp1!=itmp2) stop("Mixed type models not allowed")
	itmp<-itmp1
	tab<-merge(tab1,tab2,by=0)
	if (nrow(tab)==0) {
		warning("no common items")
		out<-list(ni=0)
	}
	if (nrow(tab)>0) {
		comuni<-tab$Row.names
		taball<-merge(tab1,tab2,by=0,all=T)
		niall<-nrow(taball)/itmp
		ni<-nrow(tab)/itmp
		if (ni==0) comuni<-"None"
		tabDff<-tab[substr(tab$Row.names,1,6)=="Dffclt",]
		tabDsc<-tab[substr(tab$Row.names,1,6)=="Dscrmn",]
		tabGss<-tab[substr(tab$Row.names,1,6)=="Gussng",]
		if (nrow(tabDff)==0 & ni>0) warning("missing difficulty parameters")
		if (nrow(tabGss)!=nrow(tabDff) & itmp==3) warning("missing guessing parameters but argument itmp=3")
		if (nrow(tabDsc)!=nrow(tabDff) & itmp==2) warning("missing discrimination parameters but argument itmp=2")
		ifelse (itmp>1,a1<-tabDsc$value1,a1<-rep(1,ni))
		ifelse (itmp>1,a2<-tabDsc$value2,a2<-rep(1,ni))
		b1<-tabDff$value1
		b2<-tabDff$value2
		ifelse (itmp==3,c1<-tabGss$value1,c1<-rep(0,ni))
		ifelse (itmp==3,c2<-tabGss$value2,c2<-rep(0,ni))

		if (method=="mean-sigma" & itmp>1) A<-sd(b2)/sd(b1)
		if (method=="mean-sigma" & itmp==1) A<-1
		if (method=="mean-mean")  A<-mean(a1)/mean(a2)
		if (method=="mean-gmean") A<-exp(sum(log(a1/a2)))^(1/ni)
		if (method=="mean-sigma" | method=="mean-mean" | method=="mean-gmean") 
			B<-mean(b2)-A*mean(b1)
		if (method=="Haebara" | method=="Stocking-Lord") {
			if (quadrature) {
				  gq<-gauss.quad.prob(nq,dist="normal")
				  ab<-gq$nodes
				  wt<-gq$weights
			}
			else {
				ab<-seq(-4,4,l=40)
				wt<-rep(1,40)
			}
			P2<-matrix(NA,length(ab),ni)
			for (i in 1:ni)
				P2[,i]<-irtp1(ab,diff=b2[i],discr=a2[i],guess=c2[i],D=D)
			par<-nlminb(start=c(1,0),objective=obj,P2=P2,ab=ab,a1=a1,b1=b1,c1=c1,met=method,itmp=itmp,wt=wt,D=D)$par

			A<-par[1]
			B<-par[2]
			b12<-A*b1+B
			a12<-a1/A
			P1<-matrix(NA,length(ab),ni)
			for (i in 1:ni)
				P1[,i]<-irtp1(ab,diff=b12[i],discr=a12[i],guess=c1[i],D=D)
		}
		if (!is.null(var1) & !is.null(var2)) {
			if (method=="mean-sigma") {
				partialA_b2<-A*sd(b2)^(-2)*(b2-mean(b2))/ni
				partialA_b1<--A*sd(b1)^(-2)*(b1-mean(b1))/ni
				if (itmp==1) {partialA_b2<-rep(0,ni)
					partialA_b1<-rep(0,ni)}
				partialB_b2<-1/ni-partialA_b2*mean(b1)
				partialB_b1<- -partialA_b1*mean(b1)-A/ni
				partialA_a2<-rep(0,ni)
				partialA_a1<-rep(0,ni)
				partialB_a2<-rep(0,ni)
				partialB_a1<-rep(0,ni)
				partialA_c2<-rep(0,ni)
				partialA_c1<-rep(0,ni)
				partialB_c2<-rep(0,ni)
				partialB_c1<-rep(0,ni)
			}
			if (method=="mean-mean") {
				partialA_b2<-rep(0,ni)
				partialA_b1<-rep(0,ni)
				partialB_b2<-rep(1/ni,ni)
				partialB_b1<-rep(-A/ni,ni)
				partialA_a2<-rep(-sum(a1)/(sum(a2)^2),ni)
				partialA_a1<-rep(1/sum(a2),ni)
				partialB_a2<--partialA_a2*mean(b1)
				partialB_a1<--partialA_a1*mean(b1)
				partialA_c2<-rep(0,ni)
				partialA_c1<-rep(0,ni)
				partialB_c2<-rep(0,ni)
				partialB_c1<-rep(0,ni)
			}
			if (method=="mean-gmean") {
				partialA_b2<-rep(0,ni)
				partialA_b1<-rep(0,ni)
				partialB_b2<-rep(1/ni,ni)
				partialB_b1<-rep(-A/ni,ni)
				partialA_a2<--1/ni*A/a2
				partialA_a1<-1/ni*A/a1
				partialB_a2<--partialA_a2*mean(b1)
				partialB_a1<--partialA_a1*mean(b1)
				partialA_c2<-rep(0,ni)
				partialA_c1<-rep(0,ni)
				partialB_c2<-rep(0,ni)
				partialB_c1<-rep(0,ni)
			}
			if (method=="Haebara") {
				P1<-t(P1)
				P2<-t(P2)
				tmp<-((c1+1)*P2+c1-2*(P2+c1+1)*P1+3*P1^2)*(P1-c1)/(1-c1)^2*(1-P1)*a1^2*D
				tmp<-apply(tmp,2,sum)*wt
				partialSIR_AB<-matrix(0,2,2)
				for (i in 1:length(ab)) 
					partialSIR_AB<- partialSIR_AB-tmp[i]*c(ab[i],1)%*%t(c((ab[i]-B)/A^2,1/A))
				if (itmp==1) {partialSIR_AB[1,]<-0
					partialSIR_AB[,1]<-0}
				abMAT<-matx(ab,ni)
				wtMAT<-matx(wt,ni)
				tmp<-(P1-c1)/(1-c1)*(1-P1)*a1*(1-P2)/(1-c2)*wtMAT
				tmp_a2<-tmp*(P2-c2)*D*(abMAT-b2)
				partialSIR_a2<-rbind(rowSums(tmp_a2*abMAT),rowSums(tmp_a2))
				tmp_b2<-tmp*(P2-c2)*D*(-a2)
				partialSIR_b2<-rbind(rowSums(tmp_b2*abMAT),rowSums(tmp_b2))
				partialSIR_c2<-rbind(rowSums(tmp*abMAT),rowSums(tmp))
				tmp_a1<-(((c1+1)*P2+c1-2*(P2+c1+1)*P1+3*P1^2)*a1/(1-c1)*
					D*(abMAT-A*b1-B)/A+P2-P1)*(P1-c1)/(1-c1)*(1-P1)*wtMAT
				partialSIR_a1<-rbind(rowSums(tmp_a1*abMAT),rowSums(tmp_a1))
				tmp_b1<-((c1+1)*P2+c1-2*(P2+c1+1)*P1+3*P1^2)*
					(P1-c1)/(1-c1)^2*(1-P1)*(-D*a1^2)*wtMAT
				partialSIR_b1<-rbind(rowSums(tmp_b1*abMAT),rowSums(tmp_b1))
				tmp_c1<-((c1+1)*P2+c1-2*(P2+c1+1)*P1+3*P1^2-(P2-P1)*(1-P1))*
					(1-P1)*a1/(1-c1)^2*wtMAT
				partialSIR_c1<-rbind(rowSums(tmp_c1*abMAT),rowSums(tmp_c1))
			}
			if (method=="Stocking-Lord") {
				P1<-t(P1)
				P2<-t(P2)
				tmp1<--(P1-c1)/(1-c1)*(1-P1)*a1
				tmp2<-P2-P1
				tmp3<-(1+c1-2*P1)/(1-c1)*a1
				tmp4<-(P1-c1)/(1-c1)*(1-P1)*a1*D
				tmp<-colSums(tmp1)*colSums(tmp4)+colSums(tmp2)*colSums(tmp3*tmp4)
				tmp<-tmp*wt
				partialSIR_AB<-matrix(0,2,2)
				for (i in 1:length(ab)) 
					partialSIR_AB<- partialSIR_AB-tmp[i]*c(ab[i],1)%*%t(c((ab[i]-B)/A^2,1/A))
				if (itmp==1) {partialSIR_AB[1,]<-0
					partialSIR_AB[,1]<-0}
				abMAT<-matrix(rep(ab,each=ni),nrow=ni)
				abMAT<-matx(ab,ni)
				wtMAT<-matx(wt,ni)
				tmp<-(P1-c1)/(1-c1)*(1-P1)*a1
				tmp<-matx(colSums(tmp),ni)*(1-P2)/(1-c2)*wtMAT
				tmp_a2<-tmp*(P2-c2)*D*(abMAT-b2)
				partialSIR_a2<-rbind(rowSums(tmp_a2*abMAT),rowSums(tmp_a2))
				tmp_b2<-tmp*(P2-c2)*D*(-a2)
				partialSIR_b2<-rbind(rowSums(tmp_b2*abMAT),rowSums(tmp_b2))
				partialSIR_c2<-rbind(rowSums(tmp*abMAT),rowSums(tmp))
				tmp1<--(P1-c1)/(1-c1)*(1-P1)*a1
				tmp2<-P2-P1
				tmp3<-(1+c1-2*P1)/(1-c1)*a1
				tmp4<-D*(abMAT-A*b1-B)/A
				tmp6<-(P1-c1)/(1-c1)*(1-P1)
				tmp_a1<-matx(colSums(tmp1),ni)*(tmp4*tmp6)+matx(colSums(tmp2),ni)*(tmp3*tmp4*tmp6)+matx(colSums(tmp2),ni)*(tmp6)
				tmp_a1<-tmp_a1*wtMAT
				partialSIR_a1<-rbind(rowSums(tmp_a1*abMAT),rowSums(tmp_a1))
				tmp61<- -tmp6*D*a1
				tmp_b1<-matx(colSums(tmp1),ni)*tmp61+matx(colSums(tmp2),ni)*(tmp3*tmp61)
				tmp_b1<-tmp_b1*wtMAT
				partialSIR_b1<-rbind(rowSums(tmp_b1*abMAT),rowSums(tmp_b1))
				tmp7<-(P1-c1)/(1-c1)*a1
				tmp8<-(1-P1)/(1-c1)
				tmp_c1<-matx(colSums(tmp1),ni)*tmp8-matx(colSums(tmp2),ni)*(tmp7*tmp8)
				tmp_c1<-tmp_c1*wtMAT
				partialSIR_c1<-rbind(rowSums(tmp_c1*abMAT),rowSums(tmp_c1))
			}

			if (method=="mean-mean" | method=="mean-sigma" | method=="mean-gmean") { 
				if (itmp==1) mat<-cbind(
					c(partialA_b1,partialA_b2),
					c(partialB_b1,partialB_b2))
				if (itmp==2) mat<-cbind(
					c(partialA_b1,partialA_a1,partialA_b2,partialA_a2),
					c(partialB_b1,partialB_a1,partialB_b2,partialB_a2))
				if (itmp==3) mat<-cbind(
					c(partialA_b1,partialA_a1,partialA_c1,partialA_b2,partialA_a2,partialA_c2),
					c(partialB_b1,partialB_a1,partialB_c1,partialB_b2,partialB_a2,partialB_c2))
			}
			if (method=="Haebara" | method=="Stocking-Lord") {
				if (itmp==1) partialSIR_gamma<-cbind(partialSIR_b1,partialSIR_b2)
				if (itmp==1) partialSIR_gamma[1,]<-0
				if (itmp==2) partialSIR_gamma<-cbind(partialSIR_b1,partialSIR_a1,partialSIR_b2,partialSIR_a2)
				if (itmp==3) partialSIR_gamma<-cbind(partialSIR_b1,partialSIR_a1,partialSIR_c1,partialSIR_b2,partialSIR_a2,partialSIR_c2)
				if (itmp==1) {invpartialSIR_AB<-partialSIR_AB
					invpartialSIR_AB[2,2]<-1/invpartialSIR_AB[2,2]}
				if (itmp==1) mat1<--invpartialSIR_AB%*%partialSIR_gamma
				if (itmp>1) mat1<--solve(partialSIR_AB)%*%partialSIR_gamma
				mat<-t(mat1)
			}
			if (ni==0) mat<-matrix(NA,2,2)
			colnames(mat)<-c("A","B")
			rownames(mat)<-c(paste(comuni,suff1,sep=""),paste(comuni,suff2,sep=""))
			if (ni>0) {
				var1<-var1[tab$Row.names,tab$Row.names]
				var2<-var2[tab$Row.names,tab$Row.names]
				rownames(var1)<-colnames(var1)<-paste(rownames(var1),suff1,sep="")
				rownames(var2)<-colnames(var2)<-paste(rownames(var2),suff2,sep="")
				var12<-blockdiag(var1,var2)
				varAB<-t(mat)%*%var12%*%mat
			}
			else { 
				varAB<-matrix(NA,2,2)
				var12<-NULL 
			}
		}
		if (is.null(var1) | is.null(var2)) {
			var12<-NULL
			mat<-NULL
			varAB<-matrix(NA,2,2)
		}
		taball$value12<-NA
		taball$value12[1:niall]<-A*taball$value1[1:niall]+B
		if (itmp>1) taball$value12[(niall+1):(2*niall)]<-taball$value1[(niall+1):(2*niall)]/A
		if (itmp==3) taball$value12[(2*niall+1):(3*niall)]<-taball$value1[(2*niall+1):(3*niall)]
		out<-list(tab1=tab1,tab2=tab2,tab=taball,var12=var12,
		partial=mat,A=A,B=B,varAB=varAB,commonitem=list(comuni),ni=ni)
	}
	out$forms<-forms
	out$method<-method
	out$itmp<-itmp
	class(out) <- "eqc"
	return(out)
}



matx<-function(vect,n) {
	rep(1,n)%x%t(vect)
}

summary.eqc <- function(object, ...)
{
	if (sum(object$ni>0)==length(object$ni)) {
		ct<-cbind(Estimate=c(object$A,object$B),StdErr=c(sqrt(diag(object$varAB))))
		rownames(ct)<-c("A","B")
	}
	else ct<-NULL
	out<-list(forms=object$forms,method=object$method,coefficients=ct)
	class(out)<-"summary.eqc"
	return(out)
}

print.summary.eqc <- function(x, ...)
{
	cat("Forms: ")
	cat(x$forms,"\n")
	cat("Method: ")
	cat(x$method,"\n")
	if (!is.null(x$coefficients)) { cat("Equating coefficients:\n")
		print(x$coefficients,digits=5)}
	else cat("no common items\n")
}


alldirec<-function(mods,method="mean-mean",all=FALSE,quadrature=TRUE,nq=30)
{
	options(warn=-1)
	nt<-length(mods)
	direclist<-list()
	k<-1
	for (i in 1:nt) {
		for (j in 1:nt) {
			if (i!=j) {
			tmp<-direc(mods[i],mods[j],suff1=paste(".",i,sep=""),suff2=paste(".",j,sep=""),method=method,quadrature=quadrature,nq=nq)
			if (tmp$ni>0 | all) {
				direclist[[k]]<-tmp
				names(direclist)[[k]]<-tmp$forms
				k<-k+1
				}
	}}}
	class(direclist)<-"eqclist"
	options(warn=0)
	return(direclist)
}


summary.eqclist <- function(object, ...)
{
	out<-list()
	for (i in 1:length(object))
		out[[i]]<-summary(object[[i]])
	class(out)<-"summary.eqclist"
	return(out)
}

print.summary.eqclist<-function(x, ...)
{
	for (i in 1:length(x)) {
		print(x[[i]])
		cat("\n\n")
	}
}


chainec<-function(r,direclist,f1=NULL,f2=NULL,pths=NULL)
{
	if (r<3) warning("r should be at least 3")
	if (is.null(pths)) {
		sel<-sapply(direclist,FUN= function(x)(x$ni!=0))
		nl<-names(direclist)[sel]
		nll<-strsplit(nl,split=".",fixed=TRUE)
		l<-data.frame(f1=sapply(nll,FUN=function(x) x[1]),f2=sapply(nll,FUN=function(x) x[2]))
		if (is.null(f1)) pths<-l
		if (!is.null(f1)) pths<-l[l$f1==f1,]
		colnames(pths)<-paste(colnames(pths),1,sep=".")
		if (r>3) {
			for (k in 1:(r-3)) {
			pths<-merge(pths,l,by.x=k+1,by.y=1)
			colnames(pths)<-paste(colnames(pths),k+1,sep=".")
			pths<-pths[,c(2:(k+1),1,k+2)]
			pths<-pths[pths[,k]!=pths[,k+2],]
			}
		}
	if (is.null(f2)) pths<-merge(pths,l,by.x=r-1,by.y=1)
	if (!is.null(f2)) pths<-merge(pths,l[l$f2==f2,],by.x=r-1,by.y=1)
	pths<-pths[,c(2:(r-1),1,r)]
	pths<-pths[pths[,r-2]!=pths[,r],] }
   
	pths<-pths[pths[,1]!=pths[,r],]
	nomi<-pths[,1]
	for (k in 2:r) nomi<-paste(nomi,pths[,k],sep=".")
	out<- vector("list", nrow(pths))
	for (j in 1:nrow(pths)) {
		ni<-c()
		A<-1
		B<-0
		partialA<-c()
		partialB<-c()
		varAll<-matrix(0,0,0)
		comuni<-list()
		missing<-FALSE
		varNULL<-FALSE
		for (k in 1:(r-1)) {
			nome<-paste(pths[j,k],pths[j,k+1],sep=".")
			link<-direclist[[nome]]
			if (k==1)   tab1<-link$tab1
			if (k==r-1) tab2<-link$tab2
			if (!is.null(link)) {
				ni<-c(ni,link$ni)
				if (link$ni!=0) {
					partialAk<-A*link$partial[,1] 
					partialA<-partialA*link$A
					partialA<-c(partialA,partialAk)

					partialBk<-B*link$partial[,1]+link$partial[,2] 
					partialB<-partialB*link$A
					partialB<-c(partialB,partialBk)

					A<-link$A*A  #A1...k=Ak-1k*A1...k-1
					B<-link$B+link$A*B #B1...k=Bkk-1+Akk-1*B1...k-1
					
					if (!is.null(link$var12)) varAll<-blockdiag(varAll,link$var12)
					if (is.null(link$var12)) varNULL<-TRUE
					comuni[[k]]<-link$commonitem[[1]]
				}
				else {
					warning("forms ",nome," have no common items\n")
					missing<-TRUE
				}
			}
			if (is.null(link)) {
				warning("link of forms ",nome," is missing\n")
				missing<-TRUE
				ni<-c(ni,0)
			}
		}

		if (!missing) {
			partialA<-tapply(partialA,names(partialA),sum)
			partialB<-tapply(partialB,names(partialB),sum)

			mat<-merge(partialA,partialB,by=0)
			nom<-mat$Row.names
			mat<-mat[,-1]
			colnames(mat)<-c("A","B")
			rownames(mat)<-nom
			mat<-as.matrix(mat)
			sel<-unique(rownames(varAll))
			varAB<-t(mat[sel,])%*%varAll[sel,sel]%*%mat[sel,]
			if(varNULL) varAB<-matrix(NA,2,2)
			#if (is.null(link$var12))
			taball<-merge(tab1,tab2,by=0,suffixes=c(.1,.2),all=T)
			niall<-nrow(taball)/link$itmp
			taball$value12<-NA
			taball$value12[1:niall]<-A*taball$value1[1:niall]+B
			if (link$itmp>1) taball$value12[(niall+1):(2*niall)]<-taball$value1[(niall+1):(2*niall)]/A
			if (link$itmp==3) taball$value12[(2*niall+1):(3*niall)]<-taball$value1[(2*niall+1):(3*niall)]
			out[[j]]$tab1<-tab1
			out[[j]]$tab2<-tab2
			out[[j]]$tab<-taball
			out[[j]]$varAll<-varAll
			out[[j]]$partial<-mat[sel,]
			out[[j]]$A<-A
			out[[j]]$B<-B
			out[[j]]$varAB<-varAB
			out[[j]]$commonitem<-comuni
		}
		out[[j]]$ni<-ni
		out[[j]]$forms<-nomi[j]
		if (!is.null(link)) out[[j]]$method<-link$method
		if (!is.null(link)) out[[j]]$itmp<-link$itmp
		else out[[j]]$method<-""
		class(out[[j]])<-"eqc"
	} 
	names(out)<-nomi
	class(out) <- "eqclist"
	return(out)
}



bisectorec<-function(ecall,mods,weighted=TRUE,unweighted=TRUE)
{
	if (length(table(sapply(ecall,FUN=function(x) x$method)))!=1) warning("ecall contains different methods")
	varNULL<-FALSE
	if (any(sapply(ecall,FUN=function(x) is.na(x$varAB)))) varNULL<-TRUE
	if (varNULL & weighted) {
		warning("weighted bisector unfeasible with NULL covariance matrix")
		weighted<-FALSE
	}
	if (!varNULL) {
		part<-lapply(ecall,FUN=function(x) data.frame(A=x$partial[,1],
					B=x$partial[,2],stringsAsFactors = FALSE))
		for (i in 1:length(part)) {
			part[[i]]$path<-names(part)[i]
			part[[i]]$par<-rownames(part[[i]])
		}
		partall<-part[[1]]
		for (i in 2:length(part))  partall<-rbind(partall,part[[i]])
		partall$link<-path2link(partall$path)
	}
	else partall<-NULL
	coall<-data.frame(t(sapply(ecall,FUN=function(x) x[c("A","B")])))
	coall$sdA<-sapply(ecall,FUN=function(x) x$varAB[1,1]^0.5)
	coall$sdB<-sapply(ecall,FUN=function(x) x$varAB[2,2]^0.5)
	for (i in 1:4) coall[,i]<-unlist(coall[,i])
	coall$path<-rownames(coall)
	coall$link<-path2link(coall$path)

	coall$weights<-NA
	links<-sort(unique(coall$link))
	if (!varNULL) VarAll<-VarExt(mods)
	else VarAll<-NULL
	if (unweighted) {
		coall$weights<-1
		bis<-bisco(coall,VarAll,partall)
	}
	if (weighted) {
		for (i in links) {
			colink<-coall[coall$link==i,]
			nl<-nrow(colink)
			weights<-rep(1,nl)
			partlink<-partall[partall$link==i,]
			o<-optim(par=weights,fn=VarTrasf,colink=colink,VarAll=VarAll,partlink=partlink,control=list(maxit=10000,reltol=1e-5))
			coall[coall$link==i,]$weights<-abs(o$par)
		}
	wbis<-bisco(coall,VarAll,partall)
	wbis$path<-"weighted bisector"
	}
	sel<-c("link","path","A","B","sdA","sdB","weights")
	if (unweighted) coall<-rbind(coall[,sel],bis[,sel])
	if (weighted) coall<-rbind(coall[,sel],wbis[,sel])
	coall<-coall[order(coall[,1]),]
	rownames(coall)<-NULL
	meq<-list(coef=coall,method=ecall[[1]]$method)
	class(meq)<-"meqc"
	return(meq)
}


summary.meqc <- function(object, ...)
{
	method<-object$method
	object<-object$coef
	link<-sort(unique(object$link))
	tab<-list()
	for (i in 1:length(link)) {
		objecti<-object[object$link==link[i],]
		objecti$coefA<-"A"
		objecti$coefB<-"B"
		objecti1<-objecti[,c("coefA","path","A","sdA")]
		objecti2<-objecti[,c("coefB","path","B","sdB")]
		colnames(objecti1)<-colnames(objecti2)<-c("","Path","Estimate","StdErr")
		tab[[i]]<-rbind(objecti1,objecti2)
	}
	out<-list(link=link,method=method,coefficients=tab)
	class(out)<-"summary.meqc"
	return(out)
}

print.summary.meqc <- function(x, ...)
{
	for (i in 1:length(x$link)) {
		cat("Link: ")
		cat(x$link[i],"\n")
		cat("Method: ")
		cat(x$method,"\n")
		cat("Equating coefficients:\n")
		print(x$coefficients[[i]],row.names=F,digits=5)
		cat("\n")
	}
}




bisco<-function(coall,VarAll,partall)
{
	coall$w<-BisW(coall$A)*coall$weights
	W<-tapply(coall$w,coall$link,FUN=sum)
	mA<-tapply(coall$A*coall$w,coall$link,sum)/W
	mB<-tapply(coall$B*coall$w,coall$link,sum)/W
	coall$W<-W[coall$link]
	coall$mA<-mA[coall$link]
	coall$mB<-mB[coall$link]
	out<-data.frame(link=names(mA),A=mA,B=mB,sdA=NA,sdB=NA,corAB=NA,path="bisector",weights=NA)
	if (!is.null(VarAll)) {
		coall$partAA<-((1-coall$A^2/(1+coall$A^2))*coall$w*coall$W+
						coall$mA*coall$W*coall$weights*coall$A*(1+coall$A^2)^(-1.5))/
						coall$W^2
		coall$partBA<- -(coall$A*coall$B/(1+coall$A^2)*coall$w*coall$W+
						coall$mB*coall$W*coall$weights*coall$A*(1+coall$A^2)^(-1.5))/
						coall$W^2
		coall$partBB<- coall$w/coall$W

		links<-unique(coall$link)
		for (i in links){
			coi<-coall[coall$link==i,]
			parti<-partall[partall$link==i,]
			parti$partAA<-coi[parti$path,]$partAA
			parti$partBA<-coi[parti$path,]$partBA
			parti$partBB<-coi[parti$path,]$partBB

			partialAmean<-tapply(parti$partAA*parti$A,parti$par,FUN=sum)
			partialBmean<-tapply(parti$partBA*parti$A+parti$partBB*parti$B,parti$par,FUN=sum)
			partialmean<-cbind(partialAmean,partialBmean)
			sel<-rownames(partialmean)
			varAB<-t(partialmean)%*%VarAll[sel,sel]%*%partialmean
			sdA<-varAB[1,1]^0.5
			sdB<-varAB[2,2]^0.5
			out[i,]$sdA<-sdA
			out[i,]$sdB<-sdB
			out[i,]$corAB<-varAB[1,2]/(sdA*sdB)
		}
	}
	return(out)
}
 
 
 
 
VarTrasf<-function(weights,colink,VarAll,partlink)
 {
	colink$weights<-abs(weights)
	colink$w<-BisW(colink$A)*colink$weights
	W<-sum(colink$w)
	mA<-sum(colink$A*colink$w)/W
	mB<-sum(colink$B*colink$w)/W
	colink$W<-W
	colink$mA<-mA
	colink$mB<-mB

	colink$partAA<-((1-colink$A^2/(1+colink$A^2))*colink$w*colink$W+
					colink$mA*colink$W*colink$weights*colink$A*(1+colink$A^2)^(-1.5))/
					colink$W^2
	colink$partBA<- -(colink$A*colink$B/(1+colink$A^2)*colink$w*colink$W+
					colink$mB*colink$W*colink$weights*colink$A*(1+colink$A^2)^(-1.5))/
					colink$W^2
	colink$partBB<- colink$w/colink$W

	partlink$partAA<-colink[partlink$path,]$partAA
	partlink$partBA<-colink[partlink$path,]$partBA
	partlink$partBB<-colink[partlink$path,]$partBB

	partialAmean<-tapply(partlink$partAA*partlink$A,partlink$par,FUN=sum)
	partialBmean<-tapply(partlink$partBA*partlink$A+partlink$partBB*partlink$B,partlink$par,FUN=sum)
	partialmean<-cbind(partialAmean,partialBmean)
	sel<-rownames(partialmean)
	varAB<-t(partialmean)%*%VarAll[sel,sel]%*%partialmean
	out<-sum(diag(varAB))
	return(out)
}



BisW<-function(x) (1+x^2)^(-0.5)


VarExt<-function(mods)
{
	VarAll<-matrix(0,0,0)
	for (i in 1:length(mods)) {
		var1<-mods[[i]]$var
		rownames(var1)<-colnames(var1)<-paste(rownames(var1),i,sep=".")
		VarAll<-blockdiag(VarAll,var1)
	}
	return(VarAll)
}



path2link<-function(x)
{
	tt<-strsplit(x,".",fixed=TRUE)
	tmp1<-sapply(tt,FUN=function(x) x[1])
	tmp2<-sapply(tt,FUN=function(x) x[length(x)])
	return(paste(tmp1,tmp2,sep="."))
}



blockdiag<-function(m1,m2)
{
	r1<-nrow(m1)
	r2<-nrow(m2)
	c1<-ncol(m1)
	c2<-ncol(m2)
	out<-matrix(0,r1+r2,c1+c2)
	if (r1>0) out[1:r1,1:c1]<-m1
	out[(r1+1):(r1+r2),(c1+1):(c1+c2)]<-m2
	rownames(out)<-c(rownames(m1),rownames(m2))
	colnames(out)<-c(colnames(m1),colnames(m2))
	return(out)
} 


convert<-function(A,B,coef=NULL,person.par=NULL)
{
	if (!is.null(coef)) {
		itms<-names(coef)
		Dffclt<-coef[substr(itms,1,6)=="Dffclt"]
		Dscrmn<-coef[substr(itms,1,6)=="Dscrmn"]
		Gussng<-coef[substr(itms,1,6)=="Gussng"]
		Dffclt<-Dffclt*A+B
		if (length(Dscrmn)>0) Dscrmn<-Dscrmn/A
		coef1<-c(Dffclt,Dscrmn,Gussng)
	}
	else coef1<-NULL
	if (!is.null(person.par)) person.par1<-person.par*A+B
	else person.par1<-NULL
	return(list(coef=coef1,person.par=person.par1))
}


import.flexmirt<-function(fnamep,fnamev=NULL,param.number=NULL,fixed=NULL,display=TRUE,digits=2) {
	par<-read.table(fnamep,fill=TRUE)
	ngr<-sum(par$V1==0,na.rm=TRUE)
	if (ngr>1) stop("Cannot handle multiple groups")
	if(max(par$V4,na.rm=TRUE)>1) stop("Cannot handle multiple factors")
	if(max(par$V6,na.rm=TRUE)>2) stop("Cannot handle multiple response models")
	par<-par[par$V1!=0 & !is.na(par$V1),]
	if(length(unique(par$V5))>1) stop("Cannot handle mixed item types")
	if(unique(par$V5)==3) stop("Cannot handle nominal categories models")
	if(unique(par$V5)==2) itmp=2
	if(unique(par$V5)==1) itmp=3
	if(length(unique(par$V8))==1) itmp=1
	ifelse(all(par$V8==1),Rasch<-TRUE,Rasch<-FALSE)
	if (itmp==1 | itmp==2) {p<-par[,c(7,8)]
		colnames(p)<-c("c","a")}
	if (itmp==3) {p<-par[,7:9]
		colnames(p)<-c("logit-g","c","a")}
	rownames(p)<-par$V2
	p<-as.matrix(p)

	if (!is.null(fnamev)) {
		vm<-read.table(fnamev,sep=",")
		vm<-as.matrix(vm[,colSums(is.na(vm))!=nrow(vm)])
		colnames(vm)<-NULL
		n<-nrow(p)
		if (itmp==1 & !Rasch) reord<-c((n+1):(2*n),1:n)
		if (itmp==1 & Rasch) reord<-1:n
		if (itmp==2) reord<-c((n+1):(2*n),1:n)
		if (itmp==3) reord<-c((2*n+1):(3*n),(n+1):(2*n),1:n)
		if (is.null(param.number)) {
			if (Rasch) param.number<-1:n
			if (itmp==1 & !Rasch) param.number<-c(rep(n+1,n),1:n)
			if (itmp==2) param.number<-c(1:n*2,1:n*2-1)
			if (itmp==3) param.number<-c(1:n*3,1:n*3-1,1:n*3-2)
		}
		if (!is.null(fixed)) {
			nf<-length(fixed)
			nvm1<-nrow(vm)+nf
			vm1<-matrix(0,nvm1,nvm1)
			vm1[!(1:nvm1)%in%fixed,!(1:nvm1)%in%fixed]<-vm
			vm<-vm1
		}
		vm<-vm[param.number,param.number]
		if (!Rasch & n*ncol(p)!=nrow(vm)) stop("number of parameters and dimension of the covariance matrix do not match")
		if (Rasch & n!=nrow(vm)) stop("number of parameters and dimension of the covariance matrix do not match")
		vm<-vm[reord,reord]
	}
	else vm<-NULL
	if (display) {
		out<-matrix(NA,nrow(p),ncol(p)*3)
		out[,seq(2,ncol(out),by=3)]<-p
		se<-diag(vm)^0.5
		if (!is.null(vm)) out[,seq(3,ncol(out),by=3)]<-se
		if (!is.null(vm)) out[,seq(1,ncol(out),by=3)]<-param.number[reord]
		if (itmp==1) colnames(out)<-c("par.num.","c","s.e.","par.num.","a","s.e.")
		if (itmp==2) colnames(out)<-c("par.num.","c","s.e.","par.num.","a","s.e.")
		if (itmp==3) colnames(out)<-c("par.num.","logit-g","s.e.","par.num.","c","s.e.","par.num.","a","s.e.")
		if (Rasch) out[,c(4,6)]<-NA
		rownames(out)<-par$V2
		out<-round(out,digits)
		print(out)
	}
	return(list(coef=p,var=vm))
}


import.irtpro<-function(fnamep,fnamev=NULL,param.number=NULL,fixed=NULL,display=TRUE,digits=2) {
	par<-read.table(fnamep,fill=TRUE)
	ngr<-sum(par$V3==0,na.rm=TRUE)
	if (ngr>1) stop("Cannot handle multiple groups")
	if(max(par$V2,na.rm=TRUE)>1) stop("Cannot handle multiple factors")
	if(max(par$V4,na.rm=TRUE)>2) stop("Cannot handle multiple response models")
	par<-par[par$V3!=0 & !is.na(par$V3),]
	if(length(unique(par$V3))>1) stop("Cannot handle mixed item types")
	if(unique(par$V3)==2) itmp=2
	if(unique(par$V3)==1) itmp=3
	if(length(unique(par$V5))==1) itmp=1
	ifelse(all(par$V5==1),Rasch<-TRUE,Rasch<-FALSE)
	if (itmp==1 | itmp==2) {p<-par[,c(6,5)]
		colnames(p)<-c("c","a")}
	if (itmp==3) {p<-par[,c(7,6,5)]
		colnames(p)<-c("g","c","a")}
	rownames(p)<-par$V1
	p<-as.matrix(p)

	if (!is.null(fnamev)) {
		vm<-read.table(fnamev,sep=",")
		vm<-as.matrix(vm[,colSums(is.na(vm))!=nrow(vm)])
		colnames(vm)<-NULL
		n<-nrow(p)
		if (itmp==1 & !Rasch) reord<-c((n+1):(2*n),1:n)
		if (itmp==1 & Rasch) reord<-1:n
		if (itmp==2) reord<-c((n+1):(2*n),1:n)
		if (itmp==3) reord<-c((2*n+1):(3*n),(n+1):(2*n),1:n)
		if (is.null(param.number)) {
			if (Rasch) param.number<-1:n
			if (itmp==1 & !Rasch) param.number<-c(rep(n+1,n),1:n)
			if (itmp==2) param.number<-c(1:n*2,1:n*2-1)
			if (itmp==3) param.number<-c(1:n*3,1:n*3-1,1:n*3-2)
		}
		if (!is.null(fixed)) {
			nf<-length(fixed)
			nvm1<-nrow(vm)+nf
			vm1<-matrix(0,nvm1,nvm1)
			vm1[!(1:nvm1)%in%fixed,!(1:nvm1)%in%fixed]<-vm
			vm<-vm1
		}
		vm<-vm[param.number,param.number]
		if (!Rasch & n*ncol(p)!=nrow(vm)) stop("number of parameters and dimension of the covariance matrix do not match")
		if (Rasch & n!=nrow(vm)) stop("number of parameters and dimension of the covariance matrix do not match")
		vm<-vm[reord,reord]
	}
	else vm<-NULL
	if (display) {
		out<-matrix(NA,nrow(p),ncol(p)*3)
		out[,seq(2,ncol(out),by=3)]<-p
		se<-diag(vm)^0.5
		if (!is.null(vm)) out[,seq(3,ncol(out),by=3)]<-se
		if (!is.null(vm)) out[,seq(1,ncol(out),by=3)]<-param.number[reord]
		if (itmp==1) colnames(out)<-c("par.num.","c","s.e.","par.num.","a","s.e.")
		if (itmp==2) colnames(out)<-c("par.num.","c","s.e.","par.num.","a","s.e.")
		if (itmp==3) colnames(out)<-c("par.num.","g","s.e.","par.num.","c","s.e.","par.num.","a","s.e.")
		if (Rasch) out[,c(4,6)]<-NA
		rownames(out)<-par$V1
		out<-round(out,digits)
		print(out)
	}
	return(list(coef=p,var=vm))
}


import.ltm<-function(mod,display=TRUE,digits=4) {
	if (class(mod)=="grm") stop("Cannot handle multiple response models")
	if (class(mod)=="gpcm") stop("Cannot handle multiple response models")
	if (class(mod)=="ltm") if (mod$ltst$factors>1) stop("Cannot handle multiple factors")
	if (class(mod)=="ltm") if (ncol(mod$coef)>2) stop("Cannot handle not IRT models")
	p<-mod$coef
	vm<-solve(mod$hessian)
	if (!is.null(mod$constraint) & class(mod)=="rasch") {
		if (!all(p[,ncol(p)]==1)) {
			cstr<-mod$constraint
			nf<-nrow(cstr)
			fixed<-cstr[,1]
			nvm1<-nrow(vm)+nf
			vm1<-matrix(0,nvm1,nvm1)
			vm1[!(1:nvm1)%in%fixed,!(1:nvm1)%in%fixed]<-vm
			vm<-vm1
		}
	}
	if (!is.null(mod$constraint) & (class(mod)=="ltm" | class(mod)=="tpm")) {
		cstr<-mod$constraint
		nf<-nrow(cstr)
		fixed<-cstr[,1]+(cstr[,2]-1)*nrow(p)
		nvm1<-nrow(vm)+nf
		vm1<-matrix(0,nvm1,nvm1)
		vm1[!(1:nvm1)%in%fixed,!(1:nvm1)%in%fixed]<-vm
		vm<-vm1
	}
	if (length(unique(p[,ncol(p)]))==1 & !all(p[,ncol(p)]==1)) {
		nv<-nrow(vm)-1
		vm<-vm[c(1:nv,rep(nv+1,nrow(p))),c(1:nv,rep(nv+1,nrow(p)))]
	}
	if (display) {
		out<-matrix(NA,nrow(p),ncol(p)*2)
		out[,seq(1,ncol(out),by=2)]<-p
		se<-diag(vm)^0.5
		if (!is.null(vm)) out[,seq(2,ncol(out),by=2)]<-se
		if (all(p[,ncol(p)]==1)) out[,ncol(out)]<-0
		rownames(out)<-rownames(p)
		colnames(out)<-1:ncol(out)
		colnames(out)[seq(1,ncol(out),by=2)]<-colnames(p)
		colnames(out)[seq(2,ncol(out),by=2)]<-rep("s.e.",ncol(p))
		out<-round(out,digits)
		print(out)
	}
	return(list(coef=p,var=vm))
}


