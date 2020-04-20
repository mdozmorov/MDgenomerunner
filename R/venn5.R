# http://bioinfo-mite.crb.wsu.edu/Rcode/Venn2.R
###An example of how to use the function
# list1=letters[1:6];
# list2=letters[6:10];
# list3=letters[6:15];
# list4=letters[6:21];
# list5=letters[6:26];
# plot_venn5(list1, list2, list3, list4, list5);
###=================================================
  
  ### here is code for the function
  
  plot_venn5<-function(listA, listB, listC, listD, listE){
    
    all_ids=unique(c(listA, listB, listC, listD, listE));
    all_n=length(all_ids);
    
    #create counting matrix
    
    all_matrix=cbind(affy=rep(0,all_n),est=rep(0,all_n),glomchip=rep(0,all_n),sage=rep(0,all_n),stanford=rep(0,all_n));
    rownames(all_matrix)=all_ids;
    colnames(all_matrix)=c(substitute(listA), substitute(listB),
                           substitute(listC), substitute(listD), substitute(listE)); 
    all_matrix[all_ids %in% listA, 1]=1;
    all_matrix[all_ids %in% listB, 2]=1;
    all_matrix[all_ids %in% listC, 3]=1;
    all_matrix[all_ids %in% listD, 4]=1;
    all_matrix[all_ids %in% listE, 5]=1;
    
    #elements only in one data set.
    
    only1= apply(all_matrix,1,sum) == 1;
    nA=sum(all_matrix[only1, 1]);
    nB=sum(all_matrix[only1, 2]);
    nC=sum(all_matrix[only1, 3]);
    nD=sum(all_matrix[only1, 4]);
    nE=sum(all_matrix[only1, 5]);
    
    #elements overlapped in two data set.
    
    only2= apply(all_matrix,1,sum) == 2;
    nAB=sum(only2 & all_matrix[,1] & all_matrix[,2]);
    nAC=sum(only2 & all_matrix[,1] & all_matrix[,3]);
    nAD=sum(only2 & all_matrix[,1] & all_matrix[,4]);
    nAE=sum(only2 & all_matrix[,1] & all_matrix[,5]);
    nBC=sum(only2 & all_matrix[,2] & all_matrix[,3]);
    nBD=sum(only2 & all_matrix[,2] & all_matrix[,4]);
    nBE=sum(only2 & all_matrix[,2] & all_matrix[,5]);
    nCD=sum(only2 & all_matrix[,3] & all_matrix[,4]);
    nCE=sum(only2 & all_matrix[,3] & all_matrix[,5]);
    nDE=sum(only2 & all_matrix[,4] & all_matrix[,5]);
    
    #elements overlapped in three data set.
    
    only3= apply(all_matrix,1,sum) == 3;
    nABC=sum(only3 & all_matrix[,1] & all_matrix[,2] & all_matrix[,3]);
    nABD=sum(only3 & all_matrix[,1] & all_matrix[,2] & all_matrix[,4]);
    nABE=sum(only3 & all_matrix[,1] & all_matrix[,2] & all_matrix[,5]);
    nACD=sum(only3 & all_matrix[,1] & all_matrix[,3] & all_matrix[,4]);
    nACE=sum(only3 & all_matrix[,1] & all_matrix[,3] & all_matrix[,5]);
    nADE=sum(only3 & all_matrix[,1] & all_matrix[,4] & all_matrix[,5]);
    nBCD=sum(only3 & all_matrix[,2] & all_matrix[,3] & all_matrix[,4]);
    nBCE=sum(only3 & all_matrix[,2] & all_matrix[,3] & all_matrix[,5]);
    nBDE=sum(only3 & all_matrix[,2] & all_matrix[,4] & all_matrix[,5]);
    nCDE=sum(only3 & all_matrix[,3] & all_matrix[,4] & all_matrix[,5]);
    
    #elements overlapped in four data set.
    
    only4= apply(all_matrix,1,sum) == 4;
    nABCD=sum(only4 & all_matrix[,1] & all_matrix[,2] & all_matrix[,3] &
                all_matrix[,4]);
    nABCE=sum(only4 & all_matrix[,1] & all_matrix[,2] & all_matrix[,3] &
                all_matrix[,5]);
    nABDE=sum(only4 & all_matrix[,1] & all_matrix[,2] & all_matrix[,4] &
                all_matrix[,5]);
    nACDE=sum(only4 & all_matrix[,1] & all_matrix[,3] & all_matrix[,4] &
                all_matrix[,5]);
    nBCDE=sum(only4 & all_matrix[,2] & all_matrix[,3] & all_matrix[,4] &
                all_matrix[,5]);
    
    #elements overlapped in five data set.
    
    all5= apply(all_matrix,1,sum) == 5;
    nABCDE=sum(all5);
    
    #make the plot.
    
    elps=cbind(150*cos(seq(0,2*pi,len=1000)), 60*sin(seq(0,2*pi,len=1000)));
    
    relocate_elp=function(e, alpha, x, y){
      phi=(alpha/180)*pi;
      xr=e[,1]*cos(phi)+e[,2]*sin(phi);
      yr=-e[,1]*sin(phi)+e[,2]*cos(phi);
      xr=x+xr;
      yr=y+yr;
      return(cbind(xr, yr));
    }
    
    par(mar=c(1,1,1,1)); 
    plot(c(0, 400), c(0, 400), type="n", axes=F, ylab="", xlab="");
    
    polygon(relocate_elp(elps, 90,200, 250));
    polygon(relocate_elp(elps, 162,250, 220));
    polygon(relocate_elp(elps, 234,250, 150));
    polygon(relocate_elp(elps, 306,180, 125));
    polygon(relocate_elp(elps, 378,145, 200));
    
    #label the data set name.
    
    text(50, 280, colnames(all_matrix)[1]);
    text(170,400, colnames(all_matrix)[2]);
    text(350,300, colnames(all_matrix)[3]);
    text(350,20, colnames(all_matrix)[4]);
    text(55,10, colnames(all_matrix)[5]);
    
    #label the numbers
    
    text(61, 228, nA);
    text(194, 329, nB);
    text(321, 245, nC);
    text(290, 81, nD);
    text(132, 69, nE);
    
    text(146, 250, nAB, cex=0.7); 
    text(123, 188, nAC, cex=0.7); 
    text(275, 152, nAD, cex=0.7); 
    text(137, 146, nAE, cex=0.7); 
    text(243, 268, nBC, cex=0.7); 
    text(175, 267, nBD, cex=0.7); 
    text(187, 117, nBE, cex=0.7); 
    text(286, 188, nCD, cex=0.7); 
    text(267, 235, nCE, cex=0.7); 
    text(228, 105, nDE, cex=0.7); 
    
    text(148, 210, nABC,cex=0.7);
    text(159, 253, nABD,cex=0.7); 
    text(171, 141, nABE,cex=0.7); 
    text(281, 175, nACD,cex=0.7); 
    text(143, 163, nACE,cex=0.7); 
    text(252, 145, nADE,cex=0.7); 
    text(205, 255, nBCD,cex=0.7); 
    text(254, 243, nBCE,cex=0.7); 
    text(211, 118, nBDE,cex=0.7); 
    text(267, 211, nCDE,cex=0.7); 
    
    text(170, 231,nABCD,cex=0.7); 
    text(158, 169,nABCE,cex=0.7); 
    text(212, 139,nABDE,cex=0.7);
    text(263, 180,nACDE,cex=0.7); 
    text(239, 232,nBCDE,cex=0.7);
    
    text(204,190,nABCDE); 
  }