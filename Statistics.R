library('ggplot2')
#####-----------------------Read in Data for statistics-------------------------------------------------#####
Datacheck<- function(fileList){
  #Datafiles = "D:\\Justin\\TH-Zeug\\Bachelor\\Abschluss\\Praxis\\Rscript_daten35\\Daten\\APG1_pvalue.csv\nD:\\Justin\\TH-Zeug\\Bachelor\\Abschluss\\Praxis\\Rscript_daten35\\Daten\\Mappe1.csv\nD:\\Justin\\TH-Zeug\\Bachelor\\Abschluss\\Praxis\\Rscript_daten35\\Daten\\Mappe2.csv\nD:\\Justin\\TH-Zeug\\Bachelor\\Abschluss\\Praxis\\Rscript_daten35\\Daten\\Mappe3.csv"
  #pth="C:\\Users\\Justin\\Desktop\\"
  
  Datafiles = fileList;
  
  i=1
  files=list();
  while(!is.na(Datafiles)){
    files[[i]]<-regmatches(Datafiles, regexpr("/n", Datafiles), invert = TRUE)[[1]][1];
    Datafiles=regmatches(Datafiles, regexpr("/n", Datafiles), invert = TRUE)[[1]][2];
    i=i+1;
  }
  i=i-1;
  pth = substring(files[[1]],1,sapply(gregexpr("/", files[[1]]), tail, 1));
  
  #write(unlist(files),file="C:\\Users\\Justin\\Desktop\\test.txt")
  
  tables=list()
  #test=1
  for(j in 2:i-1){
    tables[[j]]<-read.table(files[[j]][1], header=TRUE, sep=";",stringsAsFactors = FALSE);
    tables[[j]]$µExpr.<-as.numeric(gsub(",",".",tables[[j]]$µExpr.));
    tables[[j]]$sdExpr.<-as.numeric(gsub(",",".",tables[[j]]$sdExpr.));
    rownames(tables[[j]])=paste0(tables[[j]]$Versuchsgruppen,"_",tables[[j]]$Versuchsbedingung)
	#test=test+1;
  }
 #write.csv2(tables[[2]], file ="C:\\Users\\Justin\\Desktop\\test.csv",row.names=FALSE,sep=";", dec="," );
  sets=list()
  for(m in 1:length(rownames(tables[[1]]))){
    sets[[m]]<-tables[[1]][m,]$µExpr.;
    for(k in 2:length(tables[[1]])){
      sets[[m]][k]<-tables[[k]][m,]$µExpr.;
    }
  }
  #####-----------------------Create Histograms and qqPlots for all groups-------------------------------------------------#####
  dir.create(paste0(pth,"tmp"));
  ynameList=list();
  for(groupnum in 1:length(sets)){
    ynameList[length(ynameList)+1]=rownames(tables[[1]])[groupnum];
    png(paste0(pth,"tmp\\","Histogram_",rownames(tables[[1]])[groupnum],".png"));#,width=3.25,height=3.25,units="in",res=1200
    hist(sets[[groupnum]],col="blue",main = paste0("Histogram of µExpr for group ",rownames(tables[[1]])[groupnum]),xlab="µExpr");
    dev.off();
 # }
#  for(groupnum in 1:length(sets)){
    png(paste0(pth,"tmp\\","qqPlot_",rownames(tables[[1]])[groupnum],".png"));#,width=3.25,height=3.25,units="in",res=1200
    qqnorm(sets[[groupnum]],main = paste0("qqPlot for group ",rownames(tables[[1]])[groupnum]));
    dev.off();
  }
#  comparisonNum1<-length(unique(tables[[1]]$Versuchsbedingung));
#  comparisonNum2<-length(unique(tables[[1]]$Versuchsgruppen))-1;
#  comp=1;
#  ynameList=list();
#  for(i in 1:comparisonNum1){
#    comp1=comp+1;
#    for(j in 1:comparisonNum2){
#      ynameList[length(ynameList)+1]=paste0(tables[[1]]$Versuchsgruppen[comp1],'_',tables[[1]]$Versuchsbedingung[comp]);
#      comp1=comp1+1;
#    }
#    comp=comp+comparisonNum2+1;
#  }
  data=list();
  data[[1]]=c(length(tables),length(sets),length(ynameList));
  data=append(data,tables);
  data=append(data,sets);
  data=append(data,ynameList);
  data=append(data,pth);
  return(data);
}

#####-----------------------Statistical tests for all groups-------------------------------------------------------------#####
#####-----------------------Statistical tests for all groups-------------------------------------------------------------#####
statisticTest<-function(Dev,pair,sets,tables,compareGroups,path,cols){
   
  if(compareGroups=="All"){
	names=rownames(tables[[1]]);
	compareGroups=NULL;
	for(i in 1:(length(names)-1)){
		for(j in (i+1):length(names)){
			compareGroups=c(compareGroups,paste0(names[[i]],",",names[[j]]))
		}
	}
  }
  color=cols
  pth=path;
  stdDev=Dev;
  paired=pair;
  results=list();
  xmeans=list();
  ymeans=list();
  rnames=list();
  significance=list();
  plots=list();
  compNum = length(compareGroups);
  
  
  for(i in 1:compNum){
    xname=substring(compareGroups[i],0,regexpr(",",compareGroups[i])-1);
    yname=substring(compareGroups[i],regexpr(",",compareGroups[i])+1,nchar(compareGroups[1]));
    rnames[length(rnames)+1]=paste0(xname,"_vs_",yname);
    x=sets[[which(rownames(tables[[1]])==xname)]];
    y=sets[[which(rownames(tables[[1]])==yname)]];
    # statTable=rbind(statTable,tables[[1]][yname,1:2]);
    xmeans[length(xmeans)+1]=mean(sets[[which(rownames(tables[[1]])==xname)]]);
    ymeans[length(ymeans)+1]=mean(sets[[which(rownames(tables[[1]])==yname)]]);
	
    if(stdDev == TRUE){
      if(paired == TRUE){
        results[[length(results)+1]]=t.test(x,y,paired=TRUE)$p.value;
      } else{
        results[[length(results)+1]]=t.test(x,y,paired=FALSE)$p.value;
      }
    } else{
      if(paired == TRUE){
        results[[length(results)+1]]=wilcox.test(x,y,paired=TRUE)$p.value;
      } else{
		results[[length(results)+1]]=wilcox.test(x,y,paired=FALSE)$p.value;
      }
    }
	
    exprData=data.frame(matrix(c(sets[[which(rownames(tables[[1]])==xname)]],sets[[which(rownames(tables[[1]])==yname)]]), nrow=length(sets[[which(rownames(tables[[1]])==xname)]])*2, byrow=T))
	
    nams=c(rep(xname,length(sets[[which(rownames(tables[[1]])==xname)]])));
    nams=c(nams,rep(yname,length(sets[[which(rownames(tables[[1]])==yname)]])))
    exprData=cbind(exprData,nams);
    colnames(exprData)=c("Expression","Group")
    plots[[length(plots)+1]]=exprData;
	
  }
  #write(plots[[1]],file="C:\\Users\\Justin\\Desktop\\test1.txt")
  if(color=="Colour"){
	myplots<-lapply(plots,function(df){
    p<-ggplot( data=df,aes(x=Group,y=Expression,color=Group)) +
      geom_boxplot(outlier.colour="red", outlier.shape=8,outlier.size=4) +
      geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5)+
      ggtitle("Expression levels") +
      theme(text = element_text(size=10),
            axis.text.x = element_text(angle=45,hjust=1))+
      xlab("Experiment Group") + ylab("\u00B5Expr")
	})
  }
  else{
	myplots<-lapply(plots,function(df){
    p<-ggplot( data=df,aes(x=Group,y=Expression)) +
      geom_boxplot(outlier.colour="red", outlier.shape=8,outlier.size=4) +
      geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5)+
      ggtitle("Expression levels") +
      theme(text = element_text(size=10),
            axis.text.x = element_text(angle=45,hjust=1))+
      xlab("Experiment Group") + ylab("\u00B5Expr")
	})
  }
  
  #myplots<-ggplot( data=plots[[1]],aes(x=Group,y=Expression,color=Group)) +
  #    geom_boxplot(outlier.colour="red", outlier.shape=8,outlier.size=4) +
  #    geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5)+
  #    ggtitle("Expression levels") +
  #    theme(text = element_text(size=10),
  #          axis.text.x = element_text(angle=45,hjust=1))+
  #    xlab("Experiment Group") + ylab("\u00B5Expr")
  
  for(i in 1:compNum){
    ggsave(plot=myplots[[i]],filename=paste(pth,"tmp\\",rnames[i],".png",sep=""),dpi=600)
  }
  #write.csv2(plots[[1]],file="C:\\Users\\Justin\\Desktop\\test1.csv",row.names=FALSE,sep=";", dec="," )
  for(i in 1:length(results)){
    if(results[[1]]<=0.05){
      significance[length(significance)+1]="yes";
    }
    else{
      significance[length(significance)+1]="no";
    }
  }

  statTable=data.frame(matrix(unlist(xmeans), nrow=length(xmeans), byrow=T))
  statTable=cbind(statTable,unlist(ymeans))
  statTable=cbind(statTable,unlist(results))
  statTable=cbind(statTable,unlist(significance))
  colnames(statTable)[1]="\u00B5Expr x";
  colnames(statTable)[2]="\u00B5Expr y";
  colnames(statTable)[3]="p-value";
  colnames(statTable)[4]="Significance";
  rownames(statTable)=rnames;
  write.csv2(statTable, file =paste(pth,"tmp\\Statistics.csv",sep=""),row.names=TRUE,sep=";", dec="," );
  
 return(rnames);
}

