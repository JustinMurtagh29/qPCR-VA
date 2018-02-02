library("ddCt")
#library("reshape2")
#library("rmngb")
#library("ggplot2")
#pth="C:\\Users\\e_jmurtagh\\Desktop\\"
#control="0Gy_"
#housegene= c("TBP","RRN18S")
#sheetname1= "APG1_1.txt"
#sheet1=paste(pth,sheetname1,sep="")
#sheetname2= "APG1_2.txt"
#sheet2=paste(pth,sheetname2,sep="")
#versuchsgruppenname = "Dosis"
##versuchsgruppen = c("0Gy","2Gy","5Gy","8Gy")
#versuchsgruppen1 = c("0Gy","2Gy","5Gy")
#versuchsgruppen2 = c("0Gy","5Gy","8Gy")
#gene = "APG1"
#svname = "APG1"
#durations = c("15min","2h","24h")
#durations1 = c("15min","2h","24h")
#durations2 = c("15min","2h","24h")
#iter = 3
#plottp= "Barplot"
#resolution="300"
#label = "Expression"
#titel = "MyNewPlot"
#ylimit = c(-1,4)
#fontSize = 20
#ori = 90

#read sheet into list
readLightCycler480 <- function(file,durations,versuchsgruppen){
  #durNum = makeCounter()
  raw.data <- read.table(file, skip=4, fill=TRUE, sep="\t",header=FALSE)
  
  ## just keep columns 1 to 7
  keeps <- 1:7
  raw.data = raw.data[keeps] ## subsetting columns
  header.data <- read.table(file, skip=3, nrows=1, fill=TRUE, sep="\t", as.is=TRUE)
  header = header.data[keeps]
  colnames(raw.data) = as.character(header[1,])
  
  ## take samplenames and get information for replicates, well number and samplename
  l=list()
  repCol=NULL
  wellCol=NULL
  sampleCol=NULL
  snglist=NULL
#  capture.output(raw.data$name,file="C:\\Users\\gal0d.WS-911-1136A\\Desktop\\test.txt")
  for(i in 1:nrow(raw.data)){
    nameCol = as.character(raw.data$name[i])
    
    temp=unlist(strsplit(nameCol,"_"))
    wellCol = c(wellCol,temp[1])
    samplename = paste(temp[3:(length(temp))], collapse = '_')
    temp2 = paste(temp[2:3],collapse='_')
    samplenamegene = paste(samplename,temp2,collapse ='_',sep='_')
    snglist= c(snglist,samplenamegene)
    sampleCol = c(sampleCol, samplename)
    if(exists(samplenamegene,where=l)){
      
      l[[samplenamegene]]=l[[samplenamegene]]+1
      repCol = c(repCol,paste("Rep_",l[[samplenamegene]],sep=""))
    }else{
      repCol = c(repCol,"Rep_1")
      l[[samplenamegene]]=1
    }
    
  }
  mydata = raw.data 
  
  ## add extracted data to dataframe and reorder it to get a valid
  mydata$well = wellCol
  mydata$sample = sampleCol
  mydata$replicate = repCol
  mydata$indiv_PCR_eff = as.numeric(gsub(",",".",mydata$indiv_PCR_eff))
  test <- mydata[, c("sample","Amplicon","Cq")]
  colnames(test) = c("Sample","Detector","Ct")
  test$Ct = as.numeric(gsub(",",".",test$Ct))
  
  ## Build an InputFrame from my modified Dataframe
  ilist <- list()
  for(dur in 1:iter){
    ilist[[length(ilist)+1]] <-assign(paste("IF_",durations[[dur]],sep=""),InputFrame(test[grep(durations[[dur]], test$Sample),]))
    #durNum(increment,1)
  }  

  effVal = NULL
  effErr = NULL
  genes = NULL
  
  ## Calculate mean and standard deviation
  for(i in levels(mydata$Amplicon)){
    m = mean(mydata[which(mydata$Amplicon== i),2])
    s = sd(mydata[which(mydata$Amplicon== i),2])
    effVal = c(effVal,m)
    effErr = c(effErr,s)
    genes = c(genes,i)
  }
  names(effVal) = genes
  names(effErr) = genes
  iter = iter
  
  #Generate list to return with the calculated values
  newList <- list()
  for(j in 1:iter){
    d = durations[[j]]
    newList[[d]] <- assign(d,ilist[j])
  }
  newList[["efficiency"]] <- effVal
  newList[["efficiencyError"]] <- effErr
  return(newList)
}

#prepare data for further usage
prepdata <- function(list){
  #Create List with expression values for the ddCt method
  iter1 = iter
  inlist= list()
  for(k in 1:iter1){
    varia = durations[[k]]
    inlist[[varia]] <- list[[k]][[1]]
  }
  
  #Write efficiency and efficiency error from list into variables to be used in the ddCt method
  eff = list$efficiency
  effErr = list$efficiencyError
  
  #Create a list with the calibration samples
  caliList = list(paste(control,durations[[1]],sep=""))
  if(iter > 1){
    for(cal in 2:iter){
      caliList = c(caliList,(paste(control,durations[[cal]],sep="")))
    }
  }
  
  #Calculate expression with ddCt method and write into reslist
  resList = list()
  for(m in 1:iter1){
    resList[[length(resList)+1]] <- assign(paste("result_",durations[[m]],sep=""), ddCtExpression(inlist[[m]],
                                                                                                  algorithm == "ddCtWithE",
                                                                                                  housekeepingGenes = housegene, 
                                                                                                  calibrationSample = caliList[[m]],  
                                                                                                  efficiencies = eff, 
                                                                                                  efficiencies.error = effErr))
  }
  #Separate reslist into list with results and list with Errors
  res_durList = list()
  resErrList = list()
  for(n in 1: iter1){
    res_durList[[length(res_durList)+1]]<- assign(paste('res_',durations[[n]],sep=""), level(resList[[n]]))
    resErrList[[length(resErrList)+1]]<- assign(paste('resErr_',durations[[n]],sep=""), levelErr(resList[[n]]))
  }
  sheetlist = res_durList
  sheetlist = c(sheetlist,resErrList)
  return(sheetlist)
}

#make data frame from prepdata
makeIndata <-function(sheetlist1,durations1,versuchsgruppen1){
  #read sheet 1
  it1 = length(durations1)
  
  mid = length(sheetlist1)/2
  max = length(sheetlist1)
  mid2 = mid+1
  resdurList1 = sheetlist1[1:mid]
  resErrList1 = sheetlist1[mid2:max]
  
  plotList1 = list()
  for(o in 1:it1){
    plotList1[[o]]<- assign(paste("X",o,sep=""), rev(resdurList1[[o]][gene,]))
    if(o == 1){
      cname1 = durations1[[o]]
    }
    else{
      cname1 = c(cname1, durations1[[o]]) 
    }
  }
  
  #Create a placeholder Dataframe
  tmpi1 = list(1:length(versuchsgruppen1))
  if(it1 > 1){
    for(i in 2:it1){
      tmpi1 = c(tmpi1,list(1:length(versuchsgruppen1)))
    }
  }
  plotdata1 <- data.frame(tmpi1,row.names = versuchsgruppen1) 
  colnames(plotdata1) <- cname1
  plotdata1$Category = row.names(plotdata1)
  
  #Fill the placeholder with the according data
  for(q in 1:it1){
    for(r in 1:length(versuchsgruppen1)){
      v1 = paste(versuchsgruppen1[[r]],'_',durations1[[q]],sep='')
      
      tmpi12 <-tryCatch(plotdata1[r,q] <- plotList1[[q]][[v1]],
                        error = function(e) plotdata1[r,q] = NA)
      plotdata1[r,q]=tmpi12
    }
  }
  
  #Melt plotdata --> necessary when using ggplot2
  
  plotdata1.molten <- melt(plotdata1, value.name="Count", variable.name="Variable", na.rm=TRUE)
  #if(resErrList1[[1]][gene,])
  moltensd1 = sort(resErrList1[[1]][gene,])
  if(it1 > 1){
    for(p in 2:it1){
      moltensd1 = c(moltensd1, sort(resErrList1[[p]][gene,]))
    }
  }
  plotdata1.molten$sd = moltensd1
  
  plotdata1.molten$Variable <- factor(plotdata1.molten$Variable)
  return(plotdata1.molten)
}

#combine data from two different excelsheets
combine2Sheets <-function(sheetlist1, sheetlist2,durations,durations1,durations2,versuchsgruppen,versuchsgruppen1,versuchsgruppen2){
  
  #read sheet 1
  it1 = length(durations1)
  
  mid = length(sheetlist1)/2
  max = length(sheetlist1)
  mid2 = mid+1
  resdurList1 = sheetlist1[1:mid]
  resErrList1 = sheetlist1[mid2:max]
  
  plotList1 = list()
  for(o in 1:it1){
    plotList1[[o]]<- assign(paste("X",o,sep=""), rev(resdurList1[[o]][gene,]))
    if(o == 1){
      cname1 = durations1[[o]]
    }
    else{
      cname1 = c(cname1, durations1[[o]]) 
    }
  }
  
  #Create a placeholder Dataframe
  tmpi1 = list(1:length(versuchsgruppen1))
  if(it1 > 1){
    for(i in 2:it1){
      tmpi1 = c(tmpi1,list(1:length(versuchsgruppen1)))
    }
  }
  plotdata1 <- data.frame(tmpi1,row.names = versuchsgruppen1) 
  colnames(plotdata1) <- cname1
  plotdata1$Category = row.names(plotdata1)
  
  #Fill the placeholder with the according data
  for(q in 1:it1){
    for(r in 1:length(versuchsgruppen1)){
      v1 = paste(versuchsgruppen1[[r]],'_',durations1[[q]],sep='')
      
      tmpi12 <-tryCatch(plotdata1[r,q] <- plotList1[[q]][[v1]],
                        error = function(e) plotdata1[r,q] = NA)
      plotdata1[r,q]=tmpi12
    }
  }
  
  #Melt plotdata --> necessary when using ggplot2
  
  plotdata1.molten <- melt(plotdata1, value.name="Count", variable.name="Variable", na.rm=TRUE)
  #if(resErrList1[[1]][gene,])
  moltensd1 = sort(resErrList1[[1]][gene,])
  if(it1 > 1){
    for(p in 2:it1){
      moltensd1 = c(moltensd1, sort(resErrList1[[p]][gene,]))
    }
  }
  plotdata1.molten$sd = moltensd1
  
  plotdata1.molten$Variable <- factor(plotdata1.molten$Variable)
  
  #read sheet 2
  
  it2 = length(durations2)
  
  mid3 = length(sheetlist2)/2
  max2 = length(sheetlist2)
  mid4 = mid3+1
  resdurList2 = sheetlist2[1:mid3]
  resErrList2 = sheetlist2[mid4:max]
  
  plotList2 = list()
  for(o in 1:it2){
    plotList2[[o]]<- assign(paste("X",o,sep=""), rev(resdurList2[[o]][gene,]))
    if(o == 1){
      cname2 = durations2[[o]]
    }
    else{
      cname2 = c(cname2, durations2[[o]]) 
    }
  }
  
  #Create a placeholder Dataframe
  tmpi2 = list(1:length(versuchsgruppen2))
  if(it2 > 1){
    for(i in 2:it2){
      tmpi2 = c(tmpi2,list(1:length(versuchsgruppen2)))
    }
  }
  plotdata2 <- data.frame(tmpi2,row.names = versuchsgruppen2) 
  colnames(plotdata2) <- cname2
  plotdata2$Category = row.names(plotdata2)
  
  #Fill the placeholder with the according data
  for(q in 1:it2){
    for(r in 1:length(versuchsgruppen2)){
      v2= paste(versuchsgruppen2[[r]],'_',durations2[[q]],sep='')
      tmpi22 <-tryCatch(plotdata2[r,q] <- plotList2[[q]][[v2]],
                        error = function(e) plotdata2[r,q] = NA)
      plotdata2[r,q]=tmpi22
    }
  }
  
  #Melt plotdata --> necessary when using ggplot2
  
  plotdata2.molten <- melt(plotdata2, value.name="Count", variable.name="Variable", na.rm=TRUE)
  #if(resErrList2[[1]][gene,])
  moltensd2 = sort(resErrList2[[1]][gene,])
  if(it2 > 1){
    for(p in 2:it2){
      moltensd2 = c(moltensd2, sort(resErrList2[[p]][gene,]))
    }
  }
  plotdata2.molten$sd = moltensd2
  
  plotdata2.molten$Variable <- factor(plotdata2.molten$Variable)
  
  #create data Frame, combining sheet 1 + 2
  
  framedimension = length(durations)
  tmpi = list(1:length(versuchsgruppen))
  if(framedimension > 1){
    for(i in 2:framedimension){
      tmpi = c(tmpi,list(1:length(versuchsgruppen)))
    }
  }
  
  for(o in 1:framedimension){
    if(o == 1){
      cname = durations[[o]]
    }
    else{
      cname = c(cname, durations[[o]]) 
    }
  }
  plotdata <- data.frame(tmpi,row.names = versuchsgruppen) 
  colnames(plotdata) <- cname
  plotdata$Category = row.names(plotdata)
  
  plotdata.molten <- melt(plotdata, value.name="Count", variable.name="Variable", na.rm=TRUE)
  frameheight = length(durations)*length(versuchsgruppen)
  plotdata.molten$sd = (1:frameheight)
  plotdata.molten$Variable <- factor(plotdata.molten$Variable)
  id1 = 0
  for(q in 1:iter){
    for(r in 1:length(versuchsgruppen)){
      #v2= paste(versuchsgruppen2[[r]],'_',durations2[[q]],sep='')
      tmid1=which(plotdata1.molten$Category == versuchsgruppen[[r]])
      tmid2=which(plotdata1.molten$Variable == durations[[q]])
      id = intersect(tmid1,tmid2)
      id1 = id1+1
      if(identical(id, integer(0))){
        tmid3=which(plotdata2.molten$Category == versuchsgruppen[[r]])
        tmid4=which(plotdata2.molten$Variable == durations[[q]])
        id2 = intersect(tmid3,tmid4)
        if(identical(id2, integer(0))){
        }
        else{
          plotdata.molten$Count[id1]<-plotdata2.molten$Count[id2] 
          plotdata.molten$sd[id1]<-plotdata2.molten$sd[id2]
        }
      }
      else{
        plotdata.molten$Count[id1]<-plotdata1.molten$Count[id] 
        plotdata.molten$sd[id1]<-plotdata1.molten$sd[id] 
      }
    }
  }
  return(plotdata.molten)
}

#combine data from one excel sheet with the data resulting from combine2Sheets
combineSheets <-function(sheetlist1,durations,durations1,versuchsgruppen,versuchsgruppen1,plotdata.molten){
  
  #read sheet 1
  it1 = length(durations1)
  
  mid = length(sheetlist1)/2
  max = length(sheetlist1)
  mid2 = mid+1
  resdurList1 = sheetlist1[1:mid]
  resErrList1 = sheetlist1[mid2:max]
  
  plotList1 = list()
  for(o in 1:it1){
    plotList1[[o]]<- assign(paste("X",o,sep=""), rev(resdurList1[[o]][gene,]))
    if(o == 1){
      cname1 = durations1[[o]]
    }
    else{
      cname1 = c(cname1, durations1[[o]]) 
    }
  }
  
  #Create a placeholder Dataframe
  tmpi1 = list(1:length(versuchsgruppen1))
  if(it1 > 1){
    for(i in 2:it1){
      tmpi1 = c(tmpi1,list(1:length(versuchsgruppen1)))
    }
  }
  plotdata1 <- data.frame(tmpi1,row.names = versuchsgruppen1) 
  colnames(plotdata1) <- cname1
  plotdata1$Category = row.names(plotdata1)
  
  #Fill the placeholder with the according data
  for(q in 1:it1){
    for(r in 1:length(versuchsgruppen1)){
      v1 = paste(versuchsgruppen1[[r]],'_',durations1[[q]],sep='')
      
      tmpi12 <-tryCatch(plotdata1[r,q] <- plotList1[[q]][[v1]],
                        error = function(e) plotdata1[r,q] = NA)
      plotdata1[r,q]=tmpi12
    }
  }
  
  #Melt plotdata --> necessary when using ggplot2
  
  plotdata1.molten <- melt(plotdata1, value.name="Count", variable.name="Variable", na.rm=TRUE)
  #if(resErrList1[[1]][gene,])
  moltensd1 = sort(resErrList1[[1]][gene,])
  if(it1 > 1){
    for(p in 2:it1){
      moltensd1 = c(moltensd1, sort(resErrList1[[p]][gene,]))
    }
  }
  plotdata1.molten$sd = moltensd1
  
  plotdata1.molten$Variable <- factor(plotdata1.molten$Variable)
  
  #read new sheet data in dataFrame
  id1 = 0
  for(q in 1:iter){
    for(r in 1:length(versuchsgruppen)){
      #v2= paste(versuchsgruppen2[[r]],'_',durations2[[q]],sep='')
      tmid1=which(plotdata1.molten$Category == versuchsgruppen[[r]])
      tmid2=which(plotdata1.molten$Variable == durations[[q]])
      id = intersect(tmid1,tmid2)
      id1 = id1+1
      if(identical(id, integer(0))){
        tmid3=which(plotdata2.molten$Category == versuchsgruppen[[r]])
        tmid4=which(plotdata2.molten$Variable == durations[[q]])
        id2 = intersect(tmid3,tmid4)
        if(identical(id2, integer(0))){
        }
        else{
          plotdata.molten$Count[id1]<-plotdata2.molten$Count[id2] 
          plotdata.molten$sd[id1]<-plotdata2.molten$sd[id2]
        }
      }
      else{
        plotdata.molten$Count[id1]<-plotdata1.molten$Count[id] 
        plotdata.molten$sd[id1]<-plotdata1.molten$sd[id] 
      }
    }
  }
  return(plotdata.molten)
}

#plot the data and do the statistical Tests
plotter <- function(path,versuchsgruppenname,versuchsgruppen,ylabel,svname,plottp,resolution,plotdata.molten){
  resolution = as.numeric(resolution)
  
  # plot and facet by categories
  
  #######Creates standard bar plot ############################################
  if(plottp == "Barplot"){
    ggplot( data=plotdata.molten, aes(x=Category,y=Count, fill=Category)) +
      geom_hline(yintercept = 1)+
      geom_bar(position= position_dodge(), stat="identity") +
      geom_errorbar(aes(ymin=Count-sd, ymax=Count+sd),
                    width=.2,                    # Width of the error bars
                    position=position_dodge(.9))+
      facet_wrap( "Variable",scales="free"  )+
      ggtitle(svname) +
      xlab(versuchsgruppenname) + ylab(ylabel)
  }
  
  #######Creates bar plot with diviation from mean (1)#########################
  else{
    ggplot( data=plotdata.molten, aes(x=Category,y=Count-1, fill=Category)) +
      geom_bar(position= position_dodge(), stat="identity") +
      geom_errorbar(aes(ymin=Count-1-sd, ymax=Count-1+sd),
                    width=.2,                    # Width of the error bars
                    position=position_dodge(.9))+
      facet_wrap( "Variable",scales="free"  )+
      ggtitle(svname) +
      scale_y_continuous(breaks=-3:3, labels=-3:3 + 1)+
      xlab(versuchsgruppenname) + ylab(ylabel)
  }
  
  #Save the plot to a png file
  ggsave(filename=paste(pth,svname,".png",sep=""),dpi=resolution)
  
  # statistical test (t-test)
  options(scipen = 999)
  statistics = plotdata.molten[which(plotdata.molten$Count!=1),];
  pvalue=NULL;
  significant=NULL;
  for(i in 1:nrow(statistics)){
    
    mean = statistics$Count[i]
    sd = statistics$sd[i]
    n=3
    mu=1
    
    t <- (mean - mu) / (sd/sqrt(n))
    p=2*pt(-abs(t),df=n-1)
    p=round(x=p, digits=4)
    pvalue = c(pvalue,p)
    if(p <0.05){
      significant=c(significant,"ja")
    }else{
      significant=c(significant,"nein")
    }
  }
  
  statistics[,3] = round(x=statistics[,3],digits=4)
  statistics[,4] = round(x=statistics[,4],digits=4)
  
  statistics=cbind(statistics,pvalue,significant)
  colnames(statistics) <- c("Versuchsgruppen", "Versuchsbedingung","µ Expr.","sd Expr.","p-Wert","Signifikanz")
  
  #Write the results of the statistical test to a csv file
  write.csv2(statistics, file =paste(pth,svname,"_pvalue.csv",sep=""),row.names=FALSE,sep=";", dec="," )
}
plotadjust <- function(path,versuchsgruppenname,versuchsgruppen,ylabel,plottp,resolution,plotdata.molten,titel,ylimit,fontSize,ori){
  resolution = as.numeric(resolution)
  # plot and facet by categories
  
  #######Creates standard bar plot ############################################
  if(plottp == "Barplot"){
    ggplot( data=plotdata.molten, aes(x=Category,y=Count, fill=Category)) +
      geom_hline(yintercept = 1)+
      geom_bar(position= position_dodge(), stat="identity") +
      geom_errorbar(aes(ymin=Count-sd, ymax=Count+sd),
                    width=.2,                    # Width of the error bars
                    position=position_dodge(.9))+
      facet_wrap( "Variable",scales="free"  )+
      ggtitle(titel) +
      ylim(ylimit)+
      theme(text = element_text(size=fontSize),
            axis.text.x = element_text(angle=ori,hjust=1))+
      xlab(versuchsgruppenname) + ylab(ylabel)
  }
  
  #######Creates bar plot with diviation from mean (1)#########################
  else{
    ggplot( data=plotdata.molten, aes(x=Category,y=Count-1, fill=Category)) +
      geom_bar(position= position_dodge(), stat="identity") +
      geom_errorbar(aes(ymin=Count-1-sd, ymax=Count-1+sd),
                    width=.2,                    # Width of the error bars
                    position=position_dodge(.9))+
      facet_wrap( "Variable",scales="free"  )+
      ggtitle(titel) +
      theme(text = element_text(size=fontSize),
            axis.text.x = element_text(angle=ori,hjust=1))+
      scale_y_continuous(breaks=ylimit[1]:ylimit[2], labels=ylimit[1]:ylimit[2] + 1)+
      xlab(versuchsgruppenname) + ylab(ylabel)
  }
  
  #Save the plot to a png file
  #svname =  titel
  ggsave(filename=paste(pth,svname,".png",sep=""),dpi=resolution)
}
#list1 = readLightCycler480(sheet1,durations1)
#data1 = prepdata(list1)
#list2 = readLightCycler480(sheet2,durations2) ##add iter specific for sheet
#data2 = prepdata(list2)
#indata = combine2Sheets(data1,data2,durations,durations1,durations2,versuchsgruppen,versuchsgruppen1,versuchsgruppen2)
#plotter(pth,versuchsgruppenname,versuchsgruppen,ylabel,svname,plottp,resolution,indata)
#plotadjust(pth,versuchsgruppenname,versuchsgruppen,ylabel,plottp,resolution,indata,titel,ylimit,fontSize,ori)
