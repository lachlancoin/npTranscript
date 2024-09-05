

.processTable<-function(vcf3,sampleID=""){
  chroms=unique(vcf3$chrom)
  names(chroms)=chroms
  res_all=lapply(chroms, function(chr){
     vcf4=subset(vcf3,chrom==chr)
     strands = unique(vcf4$align_strand)
     names(strands)=strands
     lapply(strands, function(strand){
       vcf5= subset(vcf4, align_strand==strand)
       endStr = unique(vcf5$endStr)
       names(endStr) = endStr
       lapply(endStr, function(end){
         vcf6= subset(vcf5, endStr==end)
         ids=unique(vcf6$id)
         names(ids) = ids
         lapply(ids, function(id){
           vcf7 = subset(vcf6, id==id)
           nonNA = !is.na(vcf7$polyA)
           c(nrow(vcf7),length(which(nonNA)), sum(vcf7$polyA[nonNA]))
           })
       })
    })
  })
#  res_all$id = as.character(res_all$id)
  res_all
}
.convJson1<-function(json1, sampleID){
  num_cols = c("count","count_nonNA",  "sum_pA")
df_all= .merge1(lapply(json1, function(j1){
   .merge1(lapply(j1, function(j2){
     t=lapply(j2, function(j3){
       id = names(j3)
      df= cbind(data.frame(t(data.frame(j3))),id)
      names(df) = c("count","count_nonNA",  "sum_pA","code")
      df
     })
     .merge1(t,num_cols=num_cols,addName="end")
   }), num_cols=num_cols, addName="strand")
 }), num_cols=num_cols, addName="chrom")
  df_all$code = as.character(df_all$code)
  #df_all$chrom = as.character(df_all$chrom)
  #df_all$strand=as.character(df_all$strand)
  id=apply(df_all[,match(c("chrom","strand","end","code"), names(df_all))],1,paste,collapse="/")
  cbind(df_all, sampleID, id)
}


isoEnv<-R6Class("isoEnv", public = list(
  mydb="S4",
  dir="character",
  samples="character",
  isoforms="character",
  dbfile="character",
  initialize=function(dir, nme,reset = FALSE){
    self$dir = dir
    self$dbfile=paste(dir,paste0(nme,".sqlite"),sep="/")
#    rds=paste(dir,paste0(nme,".rds"), sep="/")
    self$mydb= dbConnect(RSQLite::SQLite(),self$dbfile )
    output=paste(dir,paste0(nme, ".json.gz"),sep="/")
  # self$output=output
 #  self$rds = rds
   print(paste("initialise dir",dir))
    tbls=(dbListTables(self$mydb))
    print("tables") ;
   print(tbls)
    if(reset){
      dbSendStatement(self$mydb,'DROP TABLE isoforms') #- just used for resetting
      dbSendStatement(self$mydb,'DROP TABLE samples')
      dbSendStatement(self$mydb,'DROP TABLE ip')
    }
    
    if(!("ip" %in% names(tbls))){
      df2=data.frame(matrix(0,nrow=0, ncol=4))
      names(df2) = c("sessionID","ip","sampleID","date")
      df2$ip = as.character(df2$ip)
      try(dbWriteTable(self$mydb, "ip", df2,overwrite=T,append=F))
    }
    if("samples" %in% tbls){
      ab = dbGetQuery(self$mydb, 'SELECT * from samples')
      self$samples = ab$sampleID
    }else{
      nme= c("sampleID","genomic","reference", "RNA","kit","flowcell","alignment_command","alignment_version","bin0","bin1","breakThresh")
    #  df2$ip = as.character(df2$ip)
      df2 = data.frame(matrix("",nrow=0, ncol=length(nme)))
      names(df2) = nme
      try(dbWriteTable(self$mydb, "samples", df2,overwrite=T,append=F))
     self$samples = c() 
    }
    if("isoforms" %in% tbls){
    
      ab2 = dbGetQuery(self$mydb, 'SELECT * from isoforms')
      self$isoforms = ab2$id
    }else{
      self$isoforms= c()
    }
  },
  close = function(){
    dbDisconnect(self$mydb)
    closeAllConnections()
  },
  importTable=function(vcf2, sampleID,  remote,flags){
    prev= self$register(sampleID, remote,flags)
    sessionID = prev[['sessionID']]
    tbls=(dbListTables(self$mydb))
    if("isoforms" %in% tbls){
      all_d = dbGetQuery(self$mydb, 'SELECT count(id) from isoforms where sampleID = :sampleID', list(sampleID=sampleID))
      if(all_d[1,1]>0) return(list(message="already uploaded"))
    }
    json = .processTable(vcf2, sampleID)
#    return(json1)
    self$importJSON(json,sampleID,sessionID, remote, flags)
  },
  register=function(sampleID, remote, flags){ ## upload sample data and register
    info = c("genomic","reference", "type","kit","flowcell","alignment_command","alignment_version")
   # if(length(which(!(info %in% names(flags))))>0) return(list(message="need correct info", info=info))
    all_d = dbGetQuery(self$mydb, 'SELECT * from ip where sampleID = :sampleID', list(sampleID=sampleID))
    if(nrow(all_d)>0) return(list(message="already registered", sessionID=all_d$sessionID[[1]]))
    all_d = dbGetQuery(self$mydb, 'SELECT max(sessionID) from ip')
    sessionID = all_d[1,1]+1
    if(is.na(sessionID)) sessionID=1
    df2=data.frame(matrix(0,nrow=1, ncol=4))
    names(df2) = c("sessionID","sampleID","ip","date")
    df2$sessionID = sessionID
    df2$ip = remote
    df2$sampleID = sampleID
    df2$date = date()
    try(dbWriteTable(self$mydb, "ip", df2,overwrite=F,append=T))
    flags[["sampleID"]]=sampleID
    df3 = data.frame(flags[names(flags) %in% info])
    try(dbWriteTable(self$mydb, "samples", df3,overwrite=F,append=T))
    return(list(sessionID=sessionID))
  },
  extract=function(sampleID, sessionID, remote, flags){
    all_d = dbGetQuery(self$mydb, 'SELECT * from ip where sessionID = :sessionID', list(sessionID=sessionID))
    if(nrow(all_d)==0) return(list("message"="incorrect sessionID"))
    if(all_d$sampleID[[1]]!=sampleID) return(list("message"="incorrect sessionID"))
    ab2 = dbGetQuery(self$mydb, 'SELECT count,count_nonNA, sum_pA,id from  isoforms where sampleID = :sampleID order by count DESC', list(sampleID=sampleID))
    if("top" %in% names(flags)){
      ab2 = head(ab2, flags[['top']])
    }
    ab2$pA = ab2$sum_pA/ab2$count_nonNA
    ab21 = ab2[,match( c("count","pA", "count_nonNA"),names(ab2))]
    ab21 = as.list(data.frame(t(ab21)))
    ab21=lapply(ab21, function(v){
      names(v) = c("count","pA","nonNA")
      as.list(v)
    })
    names(ab21)=ab2$id
    ab21
  },
  importJSON=function(json,sampleID, sessionID, remote, flags){ 
    all_d = dbGetQuery(self$mydb, 'SELECT * from ip where sessionID = :sessionID', list(sessionID=sessionID))
    if(nrow(all_d)==0) return(list("message"="incorrect sessionID"))
    if(all_d$sampleID[[1]]!=sampleID) return(list("message"="incorrect sessionID"))
    ## could check IP here too ? 
      tbl=  .convJson1(json,sampleID)
     
      prev_count=0; prev_pA = 0; prev_num=0; prev_pA_count=0
      tbls=(dbListTables(self$mydb))
      if("isoforms" %in% tbls){
         ab2 = dbGetQuery(self$mydb, 'SELECT id, count,count_nonNA,sum_pA from isoforms where sampleID = :sampleID ', list(sampleID=sampleID))
          novel=tbl[which(!(tbl$id %in% ab2$id)),,drop=F]
          prev_count = sum(ab2$count)
          prev_pA = sum(ab2$sum_pA)
          prev_pA_count = sum(ab2$count_nonNA)
          prev_num = nrow(ab2)
      }else{
         novel=tbl
      }
      if(nrow(novel)>0){
          try(dbWriteTable(self$mydb, "isoforms", novel,overwrite=F,append=T))
      }
       dbSendStatement(self$mydb,'UPDATE isoforms SET count= count+:count WHERE id =:id  AND sampleID = :sampleID',
                     params = as.list(tbl[,match(c("id","sampleID","count"), names(tbl))]))     
     total = prev_count + sum(tbl$count)
     total_pA_count = prev_pA_count + sum(tbl$count_nonNA)
     
     total_pA = prev_pA + sum(tbl$sum_pA)
     return(list("new"=length(novel$id), prev_num = prev_num, prev_count=prev_count, 
                 avg_pa = total_pA/total_pA_count,
                 total_reads=total,
                 new_count = sum(tbl$count)))
  }
 
))



#dir = "/home/lcoin/punim1703/Lachlan/analysis2/"
