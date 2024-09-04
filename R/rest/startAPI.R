##Rscript startAPI.R '{"params_file" : "params.json","ports" : {"8001" : "1.train_indices" }}'
##Rscript startAPI.R '{"params_file" : "params.json","ports" : {"8002" : "2.train_indices" }}'
#"dataset" : "1.train_indices", 
rm(list =ls())
library(jsonlite); 

options1<-function(arg){
  if(is.null(arg)) return(NULL)
  p1 = fromJSON(arg)
  options(p1)
  invisible(p1)
}

if(FALSE){ ## debuggin
  p1 = fromJSON('{"mem" : 1e8, "params_file" : "params.json","ports" : {"8000" : "1.train_indices"}}')
  options(p1)
}
args = commandArgs(trailingOnly=TRUE)
params_file = "params.json"
params=fromJSON("params.json"); options(params)
index=1
if(length(args)>0){
 index = as.numeric(args[1])
}

.libPaths(params$lib_paths)
#library(unix);
library(DBI); library(RSQLite); #library(GenomicAlignments); 
#library(vcfR)
#rlimit_as(getOption("mem",1e8))

invisible(lapply(params$packages, function(p) try(library(p, character.only=T))))
invisible(lapply(params$packages_server, function(p) try(library(p, character.only=T))))
#plan(multisession)

#src = c("loocObj.R","trainingEnv.R","trainObj.R", "dataObj.R","stateObj.R","modelObj.R","ypredObj.R","UDVPObj.R","helper_functions.R","dlm_lib.R","depmap_functions.R")
#src="dlm_lib.R"
#invisible(try(lapply(paste(params$path_code,"2024",src,sep="/"), function(x) {print(x);source(x)})))
#src = c("helper_functions.R","dlm_lib.R")
#invisible(try(lapply(paste(getOption("path_code"),"2024",src,sep="/"), function(x) {source(x)})))


#src1=c("bamEnv.R","bam_read_helper.R","simEnv.R","treeEnv.R")
#invisible(try(lapply(paste(params$path_code,"rest",src1,sep="/"), function(x) {print(x);source(x)})))


src1=c("isoEnv.R","helper_functions.R")
invisible(try(lapply(paste(params$path_code,"rest",src1,sep="/"), function(x) {print(x);source(x)})))


#source("bamEnv.R")
#source("bam_read_helper.R")
#source("simEnv.R")
setwd(params$working_dir)
#options(bigmemory.allow.dimnames=TRUE)


#https://spartan-ood.hpc.unimelb.edu.au/rnode/spartan-bm083.hpc.unimelb.edu.au/2189/p/7879df17/echo?msg=hello
#https://spartan-ood.hpc.unimelb.edu.au/rnode/spartan-bm083.hpc.unimelb.edu.au/32238/p/dca3307e/
#url=paste("https://spartan-ood.hpc.unimelb.edu.au/rnode/spartan-bm083.hpc.unimelb.edu.au/32238",   "p/7879df17/echo?msg=hello",sep="/")
# curl --cacert downloaded.pem  -H "Content-Type: application/json" --data '{"a":4, "b":5}' https://spartan-ood.hpc.unimelb.edu.au/rnode/spartan-bm083.hpc.unimelb.edu.au/2189/p/7879df17/sum
 #curl  -H "Content-Type: application/json" --data '{"a":4, "b":5}' "http://127.0.0.1:8000/sum"

#curl --cacert downloaded.pem -X POST  https://the-url-to-access
#all_p= fromJSON(params$data_params)  #params_dlm.json
options(params); 
#options(all_p)
#allFull = .buildAll(all_p, do_depth=F)

##following does not fill in gene names so not good

ports= getOption("urls")
#ports = paste(port, names(urls))
if(is.null(ports)) stop("ports is null")
 typen = names(ports)[[index]]

#allFull = .getGenoOnly(typen)
# all_p= options1(getOption("data_params","params1.json"))  #params_dlm.json

# allFull = .buildAll(all_p, do_depth=F)
  pr(paste(params$path_code,"rest/isoAPI.R",sep="/")) %>% pr_run(port=as.numeric(ports[[index]]),  host = "0.0.0.0")
   
