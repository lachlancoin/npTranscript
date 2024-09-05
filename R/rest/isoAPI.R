typen="all"
iso = isoEnv$new(getOption("vcf_dir",""), typen)

if(FALSE){
  remote=""
#  flags =list("remove_self"=T, count_thresh=100, prop_thresh = 0.9)
  flags = list()
  f="/home/unimelb.edu.au/lcoin/Data/npTranscript/Sequins/results_20240626221126_ref/0.readToCluster.txt.gz"
  f="/home/unimelb.edu.au/lcoin/Data/npTranscript/Sequins/results_20240626221127_dorado/0.readToCluster.txt.gz"
  
  sampleID =strsplit(rev(strsplit(f,"/")[[1]])[2],"\\.")[[1]][1]
  vcf2 =read.delim(f, header=T,sep="\t")
  vcf2=subset(vcf2,source==0 )
  iso$importTable(vcf2,sampleID, "local", flags=flags)
#  dist$importVCF(vcf3, sampleID, "local",flags=flags)
  flags = list("genomic"=TRUE,"reference"="chrIS", "type"="dRNA","kit"="NA",
               "flowcell"="MinION",
               "alignment_command"="minimap2 -y -ax splice:hq -un $fa $1.fastq  | samtools view -bS > $1.bam"
                 ,"alignment_version"="2.28-r1209")
  toJSON(list(sampleID=sampleID, flags=flags))
  list(sampleID=sampleID, flags = list());
  sessionID_=iso$register(sampleID, remote,flags)
  sessionID=sessionID_$sessionID
  json = .processTable(vcf2,sampleID)
  toJSON(list( sessionID=sessionID,sampleID=sampleID, flags=list(),json=json))
  iso$importJSON(json,sampleID, sessionID, remote, flags)
}

#in1 = read_json("test.json")
#sc = sim$score(in1$ref, in1$alt)
#write_json(sc,"out.json")

#* Echo back the input
#* @param msg The message to echo
#* @get /echo
function(msg="") {
  list(msg = paste0("The message is: '", msg, "'"))
}

#' @post /upload
function(req){
  multipart <- mime::parse_multipart(req)
# return(names(multipart$upload))
  f=multipart$upload$datapath
  vcf2 =read.delim(f, header=T,sep="\t")
  flags = list()
  print(paste("running", names(f)))
  iso$importTable(vcf2,multipart$upload$name, req$REMOTE_ADDR, flags=flags)
}



#* @post /addreads
function(req, sampleID,sessionID,flags) {
   iso$importJSON(json, sampleID, sessionID, req$REMOTE_ADDR, flags) 
}

#* @post /register
function(req, sampleID,flags) {
  #alt,sessionID, from, to,start, width, sampleID
  iso$register(sampleID, req$REMOTE_ADDR,flags)
}
