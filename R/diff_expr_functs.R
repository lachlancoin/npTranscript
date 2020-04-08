#which x is significiantly more or less than expected given y
DEgenes<-function(genenames, x,y,log=T){
  size = ceiling(sum(y))
  
  #prob = y/size
  shape1 = y 
  shape2 = size -y 
  proby = y/size;
  zeros = y==0 & x== 0
  sizex = ceiling(sum(x))
  shape1x = x 
  shape2x = sizex -x 
  
  
  probx = x/sizex
  
  pbb =matrix(NA, ncol = 2, nrow = length(x))
  
  pbb[!zeros,1] = pbetabinom.ab(x[!zeros],size = sizex,shape1 = shape1[!zeros],shape2 =shape2[!zeros],log=log)
  pbb[!zeros,2] = probx[!zeros]/proby[!zeros];
  result = pbb;
  dimnames(result) = list(genenames, c("p_lt", "ratio"))
  orders =apply(result[,1,drop=F], 2, order)
  dimnames(orders)[[1]] = genenames
  probX = probx*1e6
  probY = proby*1e6
  probX1 = (x+0.5)/sum(x+.5)
  probY1 = (y+0.5)/sum(y+.5)
  ratio1 = probX1/probY1
  #ratio1_rev = probX1/probY1
#  order()
  output =  data.frame(result,probX, probY, ratio1,x,y, genenames)
  output[order(output$p_lt),]
#  output[orders[,1],,drop=F]
}

getDescr<-function(geneneames,mart){
  attr =  c('ensembl_gene_id','description')# 'go_id') #, "name_1006", "namespace_1003") #, "definition_1006")
  filt = c('ensembl_gene_id')
  #	goids = getBM(attributes =attr,   filters = filt,   values = list(ensg) ,    mart = mart) 
  
  desc =  biomaRt::getBM(attributes=attr, filters = filt, mart = mart, values = list(genenames)) 
  desc[match(genenames, desc[,1]),]

}

getGoIDs<-function(genenames, mart){
  
 # ensg = unlist(lapply(strsplit( as.character(exons$ENSG),'\\.'), .getEl,1))
  #genenames = exons_$GENENAME.1
  attr =  c('ensembl_gene_id', 'go_id') #, "name_1006", "namespace_1003") #, "definition_1006")
  filt = c('ensembl_gene_id')
  #	goids = getBM(attributes =attr,   filters = filt,   values = list(ensg) ,    mart = mart) 
  
  goids =  biomaRt::getBM(attributes=attr, filters = filt, mart = mart, values = list(genenames)) 
  goids2 = goids[goids[,2]!="",]
  
  gn = genenames[match(goids2[,1], ensg)]
  goids = cbind(goids2, gn)
  lev_all = getlev(goids$go_id)
  dimnames(lev_all)[[1]] = goids[match(lev_all[,1], goids$go_id),3]
  goObjs = list()
  goObjs[[1]] = list(goids=goids, lev_all = lev_all )
  ns = as.factor(goids$namespace)
  lev = levels(ns)
  
  for(i in 1:length(lev)){
    goids1 = goids[goids$namespace==lev[i],]
    lev_all1 = getlev(goids1$go_id)
    dimnames(lev_all1)[[1]] = goids1[match(lev_all1[,1], goids1$go_id),3]
    goObjs[[i+1]] = list(goids = goids1, lev_all = lev_all1)
  }
  lev[lev==""]  = "blank"
  names(goObjs) = c("combined", lev)
  goObjs = goObjs[names(goObjs)!="blank"]
  goObjs
}
