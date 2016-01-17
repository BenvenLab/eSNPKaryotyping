#' Plot_Zygosity_Sinle
#'
#' Plot each SNP, without any summarization
#' @param Table The LOH table containing the output of the DeletionTable function
#' @param Organism "Human" or "Mouse"
#' @export
#' @return None

Plot_Zygosity_Sinle<-function(Table,Organism){
  
  tbl=Table
  
  if(Organism == "Human"){
    centromere_pos=c(125,93.3,91,50.4,48.4,61,59.9,45.6,49,40.2,53.7,35.8,17.9,17.6,19,36.6,24,17.2,26.5,27.5,13.2,14.7,60.6,12.5)
    centromere_pos=centromere_pos*1000000
    chr_size=c(249250621,243199373,198022430,191154276,180915260,171115067,159138663,146364022,141213431,135534747,135006516,133851895,115169878,107349540,102531392,90354753,81195210,78077248,59128983,63025520,48129895,51304566,155270560,59373566)
    lb=c(seq(from = 150, to = 10, by = -20),seq(from = 10, to = 150, by = 20))
    plot(c(1:24),rep(0,24),ylim=c(-150000000,150000000),cex.axis=0.7,xlab = "Chromsome",ylab="Position (Mbp)",xaxt = "n",yaxt="n")
    axis(2, seq(from = 150000000, to = -150000000, by = -20000000), labels=lb,cex.axis=0.7,las=1)
    mx=24
  }
  
  if(Organism == "Mouse"){
    chr_size=c(195471971,182113224,160039680,156508116,151834684,149736546,145441459,129401213,124595110,130694993,122082543,120129022,120421639,124902244,104043685,98207768  ,94987271  ,90702639  ,61431566	,171031299,91744698)
    centromere_pos=rep(200000000,21)
    lb=c(seq(from = 0, to = 200, by = 20))
    plot(c(1:21),rep(200000000,21),ylim=c(0,215000000),cex.axis=0.7,xlab = "Chromsome",ylab="Position (Mbp)",xaxt = "n",yaxt="n")
    axis(2, seq(from = 200000000, to = 0, by = -20000000), labels=lb,cex.axis=0.7,las=1)
    mx=21
  }
  
  
  for (i in 1:mx){
    chr=i
    tbl2=tbl[tbl$chr.x==chr,]
    tbl2$snp=tbl2$AF1>0
    tr=tbl2[tbl2$snp==T,]
    fl=tbl2[tbl2$snp==F,]
    # b=max(max(tr$start),max(fl$start))-min(min(tr$start),min(fl$start))
    #b=b/200
    if (dim(fl)[1]>0){
      for (j in 1:dim(fl)[1]){
        x=centromere_pos[i]-fl$start[j]
        lines(c(i,i-0.4),c(x,x),col="blue",lwd=1)
      }
    }
    if(dim(tr)[1]>0){
      for (j in 1:dim(tr)[1]){
        x=centromere_pos[i]-tr$start[j]
        lines(c(i,i+0.4),c(x,x),col="red",cex=0.05,lwd=1)
      }
    }
    
    lines(c(i,i),c(centromere_pos[i],centromere_pos[i]-chr_size[i]),lwd=10,col="black")  
    lines(c(i,i),c(centromere_pos[i],centromere_pos[i]-chr_size[i]),lwd=8,col="gray50")
    lines(c(i,i),c(centromere_pos[i],centromere_pos[i]-chr_size[i]),lwd=7,col="gray53")
    lines(c(i,i),c(centromere_pos[i],centromere_pos[i]-chr_size[i]),lwd=6,col="gray59")
    lines(c(i,i),c(centromere_pos[i],centromere_pos[i]-chr_size[i]),lwd=4,col="gray75")
    lines(c(i,i),c(centromere_pos[i],centromere_pos[i]-chr_size[i]),lwd=2,col="gray85")
    lines(c(i,i),c(centromere_pos[i],centromere_pos[i]-chr_size[i]),lwd=1,col="gray90")
    if(Organism=="Human"){
      points(i,0,pch=16,col="grey13")
      if (i<23) {text(i,centromere_pos[i]+18000000,i,cex=0.8)}
      if (i==23) {text(i,centromere_pos[i]+18000000,"X",cex=0.8)}
      if (i==24) {text(i,centromere_pos[i]+18000000,"Y",cex=0.8)}}
    
    if(Organism=="Mouse"){
      points(i,200000000,pch=16,col="grey13")
      if (i<20) {text(i,212000000,i,cex=0.8)}
      if (i==20) {text(i,212000000,"X",cex=0.8)}
      if (i==21) {text(i,212000000,"Y",cex=0.8)}}
    
    
    
  }
}