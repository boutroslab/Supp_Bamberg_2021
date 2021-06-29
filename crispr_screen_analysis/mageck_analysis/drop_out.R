pdf(file='drop_out.pdf',width=4.5,height=4.5);
gstable=read.table('drop_out.gene_summary.txt',header=T)
# 
#
# parameters
# Do not modify the variables beginning with "__"

# gstablename='__GENE_SUMMARY_FILE__'
startindex=3
# outputfile='__OUTPUT_FILE__'
targetgenelist=c("ENSG00000101182","ENSG00000140534","ENSG00000132604","ENSG00000053900","ENSG00000132646","ENSG00000162521","ENSG00000136824","ENSG00000196363","ENSG00000134014","ENSG00000136997")
# samplelabel=sub('.\\w+.\\w+$','',colnames(gstable)[startindex]);
samplelabel='DMSO1,DMSO2_vs_LIB neg.'


# You need to write some codes in front of this code:
# gstable=read.table(gstablename,header=T)
# pdf(file=outputfile,width=6,height=6)


# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")

######
# function definition

plotrankedvalues<-function(val, tglist, ...){
  
  plot(val,log='y',ylim=c(max(val),min(val)),type='l',lwd=2, ...)
  if(length(tglist)>0){
    for(i in 1:length(tglist)){
      targetgene=tglist[i];
      tx=which(names(val)==targetgene);ty=val[targetgene];
      points(tx,ty,col=colors[(i %% length(colors)) ],cex=2,pch=20)
      # text(tx+50,ty,targetgene,col=colors[i])
    }
    legend('topright',tglist,pch=20,pt.cex = 2,cex=1,col=colors)
  }
}



plotrandvalues<-function(val,targetgenelist, ...){
  # choose the one with the best distance distribution
  
  mindiffvalue=0;
  randval=val;
  for(i in 1:20){
    randval0=sample(val)
    vindex=sort(which(names(randval0) %in% targetgenelist))
    if(max(vindex)>0.9*length(val)){
      # print('pass...')
      next;
    }
    mindiffind=min(diff(vindex));
    if (mindiffind > mindiffvalue){
      mindiffvalue=mindiffind;
      randval=randval0;
      # print(paste('Diff: ',mindiffvalue))
    }
  }
  plot(randval,log='y',ylim=c(max(randval),min(randval)),pch=20,col='grey', ...)
  
  if(length(targetgenelist)>0){
    for(i in 1:length(targetgenelist)){
      targetgene=targetgenelist[i];
      tx=which(names(randval)==targetgene);ty=randval[targetgene];
      points(tx,ty,col=colors[(i %% length(colors)) ],cex=2,pch=20)
      text(tx+50,ty,targetgene,col=colors[i])
    }
  }
  
}




# set.seed(1235)



pvec=gstable[,startindex]
names(pvec)=gstable[,'id']
pvec=sort(pvec);

plotrankedvalues(pvec,targetgenelist,xlab='Genes',ylab='RRA score',main=paste('Distribution of RRA scores in \\n',samplelabel))

# plotrandvalues(pvec,targetgenelist,xlab='Genes',ylab='RRA score',main=paste('Distribution of RRA scores in \\n',samplelabel))


pvec=gstable[,startindex+1]
names(pvec)=gstable[,'id']
pvec=sort(pvec);

plotrankedvalues(pvec,targetgenelist,xlab='Genes',ylab='p value',main=paste('Distribution of p values in \\n',samplelabel))

# plotrandvalues(pvec,targetgenelist,xlab='Genes',ylab='p value',main=paste('Distribution of p values in \\n',samplelabel))



# you need to write after this code:
# dev.off()






# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(159.6157336448168,31.590151499820063,63.98608347267478),c(205.65046005989873,162.4636362847889,68.18189222498133),c(90.34852847819819,18.051515142754322,19.930091573456078),c(73.99974713452423,15.795075749910032,25.174852513839255),c(155.3134227649026,153.43787871341175,119.58054944073648),c(233.1852496913496,49.641666642574386,49.30075283960188),c(115.73216266969197,40.61590907119722,15.734282821149536),c(243.94102689113512,15.795075749910032,27.27275688999253),c(93.36014609413813,40.61590907119722,16.78323500922617),c(159.6157336448168,146.66856053487888,140.55959320226918),c(100.67407458999227,24.820833321287193,31.468565642299072),c(206.08069114789015,54.15454542826296,63.98608347267478),c(80.02298236640411,18.051515142754322,76.5735097295944),c(144.12741447712568,65.43674239248442,65.03503566075142),c(73.56951604653281,15.795075749910032,38.81123095883552),c(238.77825383523808,38.35946967835294,15.734282821149536),c(114.87170049370913,4.512878785688581,4.195808752306543),c(33.127793775339335,0.0,14.6853306330729),c(92.49968391815528,83.48825753523874,53.496561591908424),c(75.29044039849849,0.0,33.56647001845234))
targetgene="ENSG00000101182"
collabel=c("LIB","DMSO1","DMSO2")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(189.30167871622479,4.512878785688581,2.0979043761532714),c(73.56951604653281,187.2844696060761,282.168138592615),c(143.69718338913427,9.025757571377161,14.6853306330729),c(111.42985178977777,2.2564393928442903,0.0),c(283.0920558983543,216.61818171305185,150.00016289495892),c(96.37176371007807,13.53863635706574,8.391617504613086),c(73.1392849585414,0.0,26.223804701915892),c(166.92966214067096,38.35946967835294,19.930091573456078),c(83.46483107033548,13.53863635706574,5.244760940383179),c(32.267331599356496,0.0,8.391617504613086),c(60.66258340679021,54.15454542826296,27.27275688999253),c(94.22060827012098,4.512878785688581,47.2028484634486),c(131.2204818373831,185.0280302132318,35.664374394605616),c(89.48806630221534,13.53863635706574,16.78323500922617),c(195.32491394810467,9.025757571377161,12.587426256919628),c(57.65096579085027,6.76931817853287,40.90913533498879),c(103.68569220593221,216.61818171305185,219.23100730801687),c(35.278949215296436,94.77045449946019,31.468565642299072),c(221.56901031558127,293.33712106975776,281.11918640453837),c(92.49968391815528,60.92386360679584,47.2028484634486))
targetgene="ENSG00000140534"
collabel=c("LIB","DMSO1","DMSO2")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(48.185881855039035,13.53863635706574,26.223804701915892),c(59.37189014281596,63.180302999640126,108.04207537189347),c(49.90680620700471,4.512878785688581,14.6853306330729),c(283.95251807433715,178.25871203469893,133.21692788573273),c(113.15077614174345,209.84886353451898,242.30795544570285),c(22.372016575553836,2.2564393928442903,6.293713128459814),c(193.60398959613897,4.512878785688581,17.832187197302808),c(183.2784434843449,103.79621207083736,95.45464911497385),c(160.04596473280822,31.590151499820063,52.447609403831784),c(123.47632225353753,282.05492410553626,217.1331029318636),c(98.52291915003516,6.76931817853287,9.44056969268972),c(32.267331599356496,20.30795453559861,24.125900325762622),c(31.406869423373657,0.0,0.0),c(56.36027252687602,146.66856053487888,23.076948137685985),c(49.90680620700471,20.30795453559861,11.538474068842993),c(61.092814494781635,0.0,4.195808752306543),c(137.6739481572544,36.103030285508645,36.71332658268225),c(91.20899065418104,0.0,3.146856564229907),c(33.55802486333076,0.0,5.244760940383179),c(12.046470463759759,6.76931817853287,19.930091573456078))
targetgene="ENSG00000132604"
collabel=c("LIB","DMSO1","DMSO2")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(53.34865491093608,0.0,13.636378444996264),c(145.41810774109996,144.41212114203458,58.7413225322916),c(63.67420102273015,15.795075749910032,7.34266531653645),c(175.10405281250792,29.333712106975774,117.4826450645832),c(138.53441033323722,27.07727271413148,14.6853306330729),c(62.3835077587559,40.61590907119722,55.59446596806169),c(49.90680620700471,27.07727271413148,98.60150567920375),c(96.37176371007807,106.05265146368164,120.6295016288131),c(101.5345367659751,139.899242356346,68.18189222498133),c(174.6738217245165,88.00113632092732,99.6504578672804),c(128.63909530943457,2.2564393928442903,46.15389627537197),c(113.58100722973487,0.0,16.78323500922617),c(75.72067148648992,92.5140151066159,113.28683631227666),c(140.2553346852029,29.333712106975774,25.174852513839255),c(200.4876870040017,78.97537874955016,145.80435414265236),c(33.55802486333076,0.0,0.0),c(82.60436889435263,108.30909085652593,65.03503566075142),c(107.98800308584642,2.2564393928442903,2.0979043761532714),c(16.34878134367396,9.025757571377161,54.54551377998506),c(19.3603989596139,0.0,11.538474068842993))
targetgene="ENSG00000053900"
collabel=c("LIB","DMSO1","DMSO2")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(83.03459998234405,2.2564393928442903,46.15389627537197),c(179.40636369242213,230.15681807011762,182.5176807253346),c(348.0569501850587,81.23181814239445,59.79027472036824),c(292.12690874617414,56.410984821107256,100.69941005535702),c(127.77863313345173,22.564393928442904,30.419613454222436),c(252.11541756297208,194.05378778460897,212.93729417955706),c(202.20861135596738,200.82310596314184,81.81827066997758),c(76.58113366247275,45.12878785688581,56.64341815613833),c(247.8131066830579,205.3359847488304,72.37770097728786),c(92.49968391815528,78.97537874955016,53.496561591908424),c(106.26707873388072,137.6428029635017,75.52455754151777),c(207.37138441186443,49.641666642574386,17.832187197302808),c(162.19712017276532,185.0280302132318,118.53159725265984),c(187.15052327626768,78.97537874955016,38.81123095883552),c(117.88331810964907,24.820833321287193,40.90913533498879),c(140.2553346852029,76.71893935670587,56.64341815613833),c(151.87157406097126,20.30795453559861,3.146856564229907),c(40.44172227119348,0.0,2.0979043761532714),c(96.37176371007807,124.10416660643597,55.59446596806169),c(163.05758234874816,40.61590907119722,29.3706612661458))
targetgene="ENSG00000132646"
collabel=c("LIB","DMSO1","DMSO2")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(309.33615226583095,155.69431810625602,124.82531038111965),c(174.24359063652508,2.2564393928442903,8.391617504613086),c(62.3835077587559,36.103030285508645,30.419613454222436),c(63.243969934738736,27.07727271413148,28.321709078069166),c(104.11592329392363,27.07727271413148,33.56647001845234),c(157.89480929285114,227.90037867727332,63.98608347267478),c(177.25520825246502,169.23295446332176,131.11902350957945),c(125.19724660550321,250.46477260571623,292.6576604733814),c(61.52304558277305,24.820833321287193,30.419613454222436),c(72.70905387054998,20.30795453559861,45.10494408729534),c(138.1041792452458,38.35946967835294,74.47560535344114),c(202.6388424439588,47.38522724973009,94.4056969268972),c(111.42985178977777,22.564393928442904,28.321709078069166),c(85.61598651029257,0.0,10.489521880766357),c(166.0691999646881,115.07840903505881,96.50360130305049),c(121.32516681358042,124.10416660643597,191.95825041802433),c(53.34865491093608,11.282196964221452,34.615422206528976),c(77.87182692644701,0.0,1.0489521880766357),c(111.42985178977777,106.05265146368164,72.37770097728786),c(74.86020931050707,4.512878785688581,7.34266531653645))
targetgene="ENSG00000162521"
collabel=c("LIB","DMSO1","DMSO2")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(91.20899065418104,29.333712106975774,49.30075283960188),c(71.84859169456713,76.71893935670587,56.64341815613833),c(126.48793986947747,49.641666642574386,53.496561591908424),c(79.59275127841269,33.84659089266435,22.02799594960935),c(236.19686730728955,58.66742421395155,44.0559918992187),c(127.77863313345173,36.103030285508645,43.00703971114206),c(94.22060827012098,415.1848482833494,265.38490358338885),c(110.56938961379493,257.2340907842491,198.25196354648415),c(87.33691086225825,29.333712106975774,91.2588403626673),c(85.61598651029257,60.92386360679584,67.13294003690469),c(130.36001966140026,51.89810603541868,69.23084441305795),c(113.15077614174345,124.10416660643597,151.04911508303553),c(135.9530238052887,51.89810603541868,36.71332658268225),c(71.84859169456713,18.051515142754322,47.2028484634486),c(80.88344454238695,36.103030285508645,35.664374394605616),c(45.174264239099095,6.76931817853287,8.391617504613086),c(50.337037294996136,13.53863635706574,18.88113938537944),c(103.68569220593221,0.0,6.293713128459814),c(94.65083935811239,42.872348464041515,36.71332658268225),c(34.4184870393136,18.051515142754322,26.223804701915892))
targetgene="ENSG00000136824"
collabel=c("LIB","DMSO1","DMSO2")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(235.76663621929814,78.97537874955016,29.3706612661458),c(153.59249841293692,67.6931817853287,137.4127366380393),c(228.88293881143542,160.20719689194462,136.36378444996265),c(103.68569220593221,18.051515142754322,93.35674473882058),c(252.54564865096353,153.43787871341175,59.79027472036824),c(422.91715949556584,153.43787871341175,184.61558510148788),c(119.17401137362333,173.74583324901036,84.96512723420749),c(72.70905387054998,11.282196964221452,19.930091573456078),c(107.98800308584642,4.512878785688581,19.930091573456078),c(156.17388494088544,65.43674239248442,54.54551377998506),c(182.41798130836207,78.97537874955016,47.2028484634486),c(135.52279271729728,4.512878785688581,8.391617504613086),c(160.47619582079963,99.28333328514877,94.4056969268972),c(93.36014609413813,24.820833321287193,22.02799594960935),c(100.67407458999227,54.15454542826296,40.90913533498879),c(162.62735126075674,257.2340907842491,140.55959320226918),c(75.72067148648992,51.89810603541868,36.71332658268225),c(95.08107044610381,151.18143932056745,96.50360130305049),c(132.51117510135734,24.820833321287193,32.51751783037571),c(63.67420102273015,18.051515142754322,48.251800651525244))
targetgene="ENSG00000196363"
collabel=c("LIB","DMSO1","DMSO2")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(264.16188802673184,33.84659089266435,58.7413225322916),c(238.77825383523808,433.2363634261037,218.18205511994023),c(160.47619582079963,83.48825753523874,58.7413225322916),c(164.77850670071385,27.07727271413148,13.636378444996264),c(183.2784434843449,273.0291665341591,244.40585982185613),c(107.98800308584642,29.333712106975774,33.56647001845234),c(80.45321345439554,24.820833321287193,58.7413225322916),c(61.953276670764474,24.820833321287193,39.860183146912156),c(167.35989322866237,103.79621207083736,56.64341815613833),c(119.17401137362333,27.07727271413148,17.832187197302808),c(140.2553346852029,119.59128782074738,50.34970502767851),c(143.26695230114285,0.0,15.734282821149536),c(51.197499470978975,4.512878785688581,5.244760940383179),c(73.99974713452423,130.87348478496884,60.83922690844487),c(65.82535646268725,63.180302999640126,32.51751783037571),c(121.75539790157185,0.0,0.0),c(195.32491394810467,90.25757571377162,57.692370344214964),c(52.91842382294465,20.30795453559861,4.195808752306543),c(67.97651190264435,6.76931817853287,2.0979043761532714),c(34.84871812730502,22.564393928442904,40.90913533498879))
targetgene="ENSG00000134014"
collabel=c("LIB","DMSO1","DMSO2")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(132.08094401336592,78.97537874955016,41.95808752306543),c(247.8131066830579,178.25871203469893,103.84626661958694),c(82.60436889435263,49.641666642574386,69.23084441305795),c(200.05745591601027,24.820833321287193,34.615422206528976),c(138.53441033323722,160.20719689194462,113.28683631227666),c(154.8831916769112,101.53977267799306,75.52455754151777),c(175.53428390049933,11.282196964221452,27.27275688999253),c(205.22022897190732,56.410984821107256,54.54551377998506),c(112.2903139657606,142.15568174919028,39.860183146912156),c(200.4876870040017,36.103030285508645,18.88113938537944),c(127.3484020454603,29.333712106975774,46.15389627537197),c(96.80199479806949,29.333712106975774,29.3706612661458),c(72.70905387054998,18.051515142754322,24.125900325762622),c(193.17375850814756,47.38522724973009,99.6504578672804),c(81.31367563037837,106.05265146368164,92.30779255074394),c(35.70918030328786,11.282196964221452,12.587426256919628),c(35.70918030328786,6.76931817853287,7.34266531653645),c(84.32529324631831,74.46249996386157,41.95808752306543),c(80.02298236640411,6.76931817853287,51.39865721575515),c(43.45333988713342,83.48825753523874,7.34266531653645))
targetgene="ENSG00000136997"
collabel=c("LIB","DMSO1","DMSO2")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}



# 
#
# parameters
# Do not modify the variables beginning with "__"

# gstablename='__GENE_SUMMARY_FILE__'
startindex=9
# outputfile='__OUTPUT_FILE__'
targetgenelist=c("ENSG00000173011","ENSG00000184634","ENSG00000183337","ENSG00000204371","ENSG00000181090","ENSG00000171148","ENSG00000152382","ENSG00000189221","ENSG00000166503","ENSG00000171843")
# samplelabel=sub('.\\w+.\\w+$','',colnames(gstable)[startindex]);
samplelabel='DMSO1,DMSO2_vs_LIB pos.'


# You need to write some codes in front of this code:
# gstable=read.table(gstablename,header=T)
# pdf(file=outputfile,width=6,height=6)


# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")

######
# function definition

plotrankedvalues<-function(val, tglist, ...){
  
  plot(val,log='y',ylim=c(max(val),min(val)),type='l',lwd=2, ...)
  if(length(tglist)>0){
    for(i in 1:length(tglist)){
      targetgene=tglist[i];
      tx=which(names(val)==targetgene);ty=val[targetgene];
      points(tx,ty,col=colors[(i %% length(colors)) ],cex=2,pch=20)
      # text(tx+50,ty,targetgene,col=colors[i])
    }
    legend('topright',tglist,pch=20,pt.cex = 2,cex=1,col=colors)
  }
}



plotrandvalues<-function(val,targetgenelist, ...){
  # choose the one with the best distance distribution
  
  mindiffvalue=0;
  randval=val;
  for(i in 1:20){
    randval0=sample(val)
    vindex=sort(which(names(randval0) %in% targetgenelist))
    if(max(vindex)>0.9*length(val)){
      # print('pass...')
      next;
    }
    mindiffind=min(diff(vindex));
    if (mindiffind > mindiffvalue){
      mindiffvalue=mindiffind;
      randval=randval0;
      # print(paste('Diff: ',mindiffvalue))
    }
  }
  plot(randval,log='y',ylim=c(max(randval),min(randval)),pch=20,col='grey', ...)
  
  if(length(targetgenelist)>0){
    for(i in 1:length(targetgenelist)){
      targetgene=targetgenelist[i];
      tx=which(names(randval)==targetgene);ty=randval[targetgene];
      points(tx,ty,col=colors[(i %% length(colors)) ],cex=2,pch=20)
      text(tx+50,ty,targetgene,col=colors[i])
    }
  }
  
}




# set.seed(1235)



pvec=gstable[,startindex]
names(pvec)=gstable[,'id']
pvec=sort(pvec);

plotrankedvalues(pvec,targetgenelist,xlab='Genes',ylab='RRA score',main=paste('Distribution of RRA scores in \\n',samplelabel))

# plotrandvalues(pvec,targetgenelist,xlab='Genes',ylab='RRA score',main=paste('Distribution of RRA scores in \\n',samplelabel))


pvec=gstable[,startindex+1]
names(pvec)=gstable[,'id']
pvec=sort(pvec);

plotrankedvalues(pvec,targetgenelist,xlab='Genes',ylab='p value',main=paste('Distribution of p values in \\n',samplelabel))

# plotrandvalues(pvec,targetgenelist,xlab='Genes',ylab='p value',main=paste('Distribution of p values in \\n',samplelabel))



# you need to write after this code:
# dev.off()






# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(125.62747769349463,239.18257564149476,262.2380470191589),c(249.10379994703214,602.4693178894255,344.0563176891365),c(114.44146940571771,428.72348464041517,388.11230958835523),c(244.37125797912654,995.089772244332,250.69957295031594),c(92.92991500614671,182.7715908203875,75.52455754151777),c(227.16201445946973,1640.4314385977991,725.874914149032),c(56.36027252687602,329.4401513552664,354.54583956990285),c(96.80199479806949,717.5477269244843,461.5389627537197),c(184.13890566032774,349.748105890865,282.168138592615),c(49.90680620700471,47.38522724973009,287.4128995329982),c(174.6738217245165,306.8757574268235,334.61574799644677),c(57.220734702858856,180.51515142754323,139.51064101419254),c(92.92991500614671,381.33825739068504,201.39882011071404),c(136.38325489328014,88.00113632092732,83.91617504613086),c(138.53441033323722,548.3147724611625,884.2666945486039),c(74.42997822251566,194.05378778460897,122.72740600496638),c(75.72067148648992,51.89810603541868,70.27979660113459),c(97.66245697405233,566.3662876039168,493.00752839601876),c(123.90655334152895,541.5454542826296,412.2382099141178),c(58.51142796683311,83.48825753523874,78.67141410574767))
targetgene="ENSG00000173011"
collabel=c("LIB","DMSO1","DMSO2")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(197.04583830007033,284.31136349838056,251.74852513839258),c(74.42997822251566,189.54090899892037,227.62262481262994),c(87.33691086225825,412.9284088905051,320.97936955145053),c(66.6858186386701,101.53977267799306,129.0211191334262),c(130.36001966140026,737.8556814600829,585.3153209467628),c(72.27882278255855,130.87348478496884,256.9932860787757),c(49.90680620700471,38.35946967835294,55.59446596806169),c(213.82485073173572,796.5231056740345,622.0286475294449),c(218.1271616116499,189.54090899892037,285.31499515684493),c(147.56926318105704,166.9765150704775,310.4898476706842),c(101.10430567798369,42.872348464041515,167.83235009226172),c(157.03434711686828,97.02689389230449,121.67845381688974),c(159.18550255682538,365.543181640775,276.9233776522318),c(55.069579262901755,117.3348484279031,115.38474068842993),c(119.60424246161475,56.410984821107256,71.32874878921123),c(157.89480929285114,749.1378784243044,286.36394734492154),c(52.48819273495324,56.410984821107256,29.3706612661458),c(178.9761326044307,555.0840906396954,402.7976402214281),c(55.49981035089318,291.08068167691346,158.39178039957199),c(38.290566831236376,133.12992417781314,279.0212820283851))
targetgene="ENSG00000184634"
collabel=c("LIB","DMSO1","DMSO2")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(133.80186836533161,81.23181814239445,182.5176807253346),c(101.10430567798369,99.28333328514877,140.55959320226918),c(183.2784434843449,324.9272725695778,271.67861671184863),c(203.0690735319502,194.05378778460897,242.30795544570285),c(166.4994310526795,539.2890148897853,439.51096680411035),c(92.06945283016387,146.66856053487888,76.5735097295944),c(117.88331810964907,117.3348484279031,284.26604296876826),c(188.44121654024195,198.56666657029754,230.76948137685986),c(210.81323311579578,681.4446966389756,403.84659240950475),c(90.34852847819819,248.20833321287193,276.9233776522318),c(171.23197302058514,746.8814390314601,352.4479351937496),c(109.7089274378121,672.4189390675986,841.2596548374619),c(108.84846526182925,78.97537874955016,111.18893193612338),c(48.185881855039035,36.103030285508645,74.47560535344114),c(82.17413780636122,383.59469678352934,227.62262481262994),c(62.3835077587559,42.872348464041515,142.65749757842246),c(94.65083935811239,11.282196964221452,29.3706612661458),c(53.77888599892749,13.53863635706574,80.76931848190095),c(109.7089274378121,424.21060585472657,238.1121466933963),c(0.0,0.0,0.0))
targetgene="ENSG00000183337"
collabel=c("LIB","DMSO1","DMSO2")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(184.56913674831915,478.36515128298953,210.83938980340378),c(234.9061740433153,968.0124995302006,300.00032578991784),c(231.89455642737536,947.704544994602,641.9587391029011),c(237.05732948327238,392.6204543549065,248.60166857416266),c(233.61548077934103,548.3147724611625,319.93041736337386),c(197.9063004760532,480.62159067583383,259.091190454929),c(33.988255951322174,200.82310596314184,195.10510698225423),c(90.7787595661896,45.12878785688581,174.12606322072153),c(94.22060827012098,121.84772721359168,106.99312318381685),c(148.85995644503132,63.180302999640126,113.28683631227666),c(94.22060827012098,88.00113632092732,67.13294003690469),c(92.92991500614671,108.30909085652593,69.23084441305795),c(229.31316989942684,735.5992420672386,420.62982741873094),c(103.2554611179408,15.795075749910032,76.5735097295944),c(92.06945283016387,27.07727271413148,101.74836224343366),c(148.85995644503132,787.4973481026573,593.7069384513758),c(67.97651190264435,203.0795453559861,181.468728537258),c(69.26720516661861,124.10416660643597,167.83235009226172),c(121.75539790157185,412.9284088905051,95.45464911497385),c(46.46495750307336,11.282196964221452,19.930091573456078))
targetgene="ENSG00000204371"
collabel=c("LIB","DMSO1","DMSO2")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(89.48806630221534,99.28333328514877,74.47560535344114),c(129.9297885734088,207.5924241416747,180.41977634918135),c(171.66220410857656,629.546590603557,319.93041736337386),c(325.2547025215135,361.03030285508646,238.1121466933963),c(141.54602794917716,216.61818171305185,158.39178039957199),c(84.75552433430973,442.2621209974809,151.04911508303553),c(156.60411602887686,410.6719694976608,239.16109888147295),c(141.97625903716857,144.41212114203458,189.86034604187105),c(136.38325489328014,478.36515128298953,309.4408954826075),c(205.22022897190732,719.8041663173286,549.6509465521572),c(114.87170049370913,365.543181640775,489.86067183178886),c(129.4995574854174,243.69545442718336,527.6229506025478),c(99.81361241400943,250.46477260571623,86.01407942228413),c(160.90642690879108,214.36174232020758,117.4826450645832),c(174.6738217245165,221.13106049874045,140.55959320226918),c(157.89480929285114,480.62159067583383,267.48280795954213),c(204.3597667959245,523.4939391398754,266.43385577146546),c(144.98787665310851,241.43901503433906,272.7275688999253),c(96.80199479806949,130.87348478496884,104.89521880766357),c(62.3835077587559,69.949621178173,93.35674473882058))
targetgene="ENSG00000181090"
collabel=c("LIB","DMSO1","DMSO2")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(95.08107044610381,245.95189382002764,104.89521880766357),c(157.4645782048597,313.64507560535634,228.67157700070658),c(71.84859169456713,189.54090899892037,75.52455754151777),c(69.69743625461004,137.6428029635017,93.35674473882058),c(128.20886422144315,528.006817925564,394.406022716815),c(166.92966214067096,162.4636362847889,180.41977634918135),c(157.03434711686828,250.46477260571623,142.65749757842246),c(155.74365385289403,254.9776513914048,238.1121466933963),c(114.44146940571771,148.92499992772315,116.43369287650657),c(274.0572030505345,1137.2454539935222,378.6717398956655),c(90.7787595661896,160.20719689194462,137.4127366380393),c(79.59275127841269,151.18143932056745,106.99312318381685),c(158.75527146883397,286.56780289122486,456.29420181333654),c(134.66233054131445,115.07840903505881,117.4826450645832),c(221.56901031558127,1031.1928025298407,355.5947917579795),c(168.2203554046452,525.7503785327197,194.05615479417762),c(91.20899065418104,254.9776513914048,136.36378444996265),c(152.30180514896267,173.74583324901036,131.11902350957945),c(129.4995574854174,69.949621178173,74.47560535344114),c(56.36027252687602,83.48825753523874,63.98608347267478))
targetgene="ENSG00000171148"
collabel=c("LIB","DMSO1","DMSO2")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(135.9530238052887,157.95075749910032,243.35690763377949),c(129.9297885734088,550.5712118540068,610.490173460602),c(163.48781343673957,728.8299238887058,612.5880778367552),c(79.59275127841269,336.2094695337993,267.48280795954213),c(158.32504038084255,108.30909085652593,147.90225851880564),c(123.90655334152895,128.61704539212454,250.69957295031594),c(206.5109222358816,616.0079542464913,286.36394734492154),c(85.18575542230116,275.2856059270034,190.9092982299477),c(96.80199479806949,139.899242356346,203.49672448686732),c(67.11604972666152,81.23181814239445,100.69941005535702),c(54.209117086918916,90.25757571377162,139.51064101419254),c(110.13915852580351,209.84886353451898,206.64358105109724),c(189.7319098042162,273.0291665341591,247.55271638608602),c(86.90667977426683,221.13106049874045,186.71348947764116),c(186.29006110028484,230.15681807011762,253.84642951454583),c(129.9297885734088,273.0291665341591,252.7974773264692),c(76.15090257448134,110.56553024937023,84.96512723420749),c(103.2554611179408,58.66742421395155,90.20988817459067),c(56.36027252687602,110.56553024937023,196.15405917033087),c(94.65083935811239,76.71893935670587,48.251800651525244))
targetgene="ENSG00000152382"
collabel=c("LIB","DMSO1","DMSO2")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(93.36014609413813,112.82196964221451,113.28683631227666),c(92.49968391815528,153.43787871341175,247.55271638608602),c(104.97638546990648,88.00113632092732,49.30075283960188),c(110.13915852580351,45.12878785688581,38.81123095883552),c(187.5807543642591,196.31022717745324,139.51064101419254),c(51.6277305589704,38.35946967835294,26.223804701915892),c(131.2204818373831,99.28333328514877,224.47576824840004),c(110.99962070178636,266.2598483556263,150.00016289495892),c(110.99962070178636,331.6965907481107,229.72052918878322),c(114.44146940571771,171.48939385616606,98.60150567920375),c(76.58113366247275,56.410984821107256,67.13294003690469),c(83.8950621583269,148.92499992772315,114.33578850035329),c(105.40661655789789,137.6428029635017,103.84626661958694),c(100.24384350200086,185.0280302132318,140.55959320226918),c(85.61598651029257,207.5924241416747,250.69957295031594),c(67.97651190264435,342.9787877123321,300.00032578991784),c(107.98800308584642,356.51742406939786,145.80435414265236),c(91.20899065418104,320.41439378388924,309.4408954826075),c(112.2903139657606,568.6227269967611,276.9233776522318),c(101.5345367659751,234.6696968558062,231.8184335649365))
targetgene="ENSG00000189221"
collabel=c("LIB","DMSO1","DMSO2")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(56.36027252687602,13.53863635706574,26.223804701915892),c(243.94102689113512,661.1367421033771,322.02832173952714),c(155.3134227649026,157.95075749910032,111.18893193612338),c(88.1973730382411,81.23181814239445,70.27979660113459),c(54.63934817491033,128.61704539212454,54.54551377998506),c(185.85983001229343,399.3897725334394,236.01424231724303),c(115.30193158170054,42.872348464041515,88.1119837984374),c(187.5807543642591,550.5712118540068,345.1052698772132),c(59.37189014281596,83.48825753523874,28.321709078069166),c(147.13903209306562,189.54090899892037,125.87426256919629),c(120.46470463759759,250.46477260571623,264.3359513953122),c(63.67420102273015,153.43787871341175,166.78339790418508),c(120.03447354960616,169.23295446332176,109.09102755997012),c(124.33678442952036,275.2856059270034,143.7064497664991),c(136.38325489328014,363.28674224793076,164.6854935280318),c(83.03459998234405,182.7715908203875,125.87426256919629),c(156.60411602887686,277.5420453198477,215.0351985557103),c(79.59275127841269,157.95075749910032,152.09806727111217),c(61.953276670764474,54.15454542826296,26.223804701915892),c(77.4415958384556,180.51515142754323,121.67845381688974))
targetgene="ENSG00000166503"
collabel=c("LIB","DMSO1","DMSO2")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(251.25495538698925,358.77386346224216,253.84642951454583),c(46.46495750307336,49.641666642574386,114.33578850035329),c(183.70867457233632,559.596969425384,223.4268160603234),c(100.67407458999227,173.74583324901036,121.67845381688974),c(67.11604972666152,115.07840903505881,203.49672448686732),c(218.1271616116499,200.82310596314184,248.60166857416266),c(71.41836060657572,291.08068167691346,174.12606322072153),c(86.04621759828399,227.90037867727332,137.4127366380393),c(123.90655334152895,232.4132574629619,166.78339790418508),c(217.26669943566708,311.38863621251204,241.2590032576262),c(97.23222588606092,162.4636362847889,112.23788412420002),c(97.23222588606092,625.0337118178684,344.0563176891365),c(67.97651190264435,60.92386360679584,36.71332658268225),c(95.94153262208665,60.92386360679584,73.4266531653645),c(147.99949426904845,194.05378778460897,124.82531038111965),c(132.94140618934875,139.899242356346,280.07023421646176),c(53.34865491093608,24.820833321287193,170.9792066564916),c(65.82535646268725,315.90151499820064,226.57367262455332),c(67.97651190264435,94.77045449946019,130.07007132150284),c(93.36014609413813,196.31022717745324,161.5386369638019))
targetgene="ENSG00000171843"
collabel=c("LIB","DMSO1","DMSO2")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}



dev.off()
Sweave("drop_out_summary.Rnw");
library(tools);

texi2dvi("drop_out_summary.tex",pdf=TRUE);

