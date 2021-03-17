package npTranscript.cluster;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.SortedMap;
import java.util.SortedSet;
import java.util.Stack;
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.ThreadPoolExecutor;

import htsjdk.samtools.SAMRecord;
import japsa.seq.Sequence;
import npTranscript.run.ViralTranscriptAnalysisCmd2;

public class IdentityProfileHolder {
	
	public static void waitOnThreads(int sleep) {
		if(executor==null) return;
		if(executor instanceof ThreadPoolExecutor){
	    	while(((ThreadPoolExecutor)executor).getActiveCount()>0){
	    		try{
		    	System.err.println("IdentityProfileHolder: awaiting completion "+((ThreadPoolExecutor)executor).getActiveCount());
		    	//Thread.currentThread();
				Thread.sleep(sleep);
	    		}catch(InterruptedException exc){
	    			exc.printStackTrace();
	    		}
	    	}
	    	}
		//executor.a
	}
	
	static int primelen = 500;//chr.length();
	public static  ExecutorService executor ;
 final CigarClusters all_clusters;
 final BreakPoints bp;
 final Outputs o;
 final Sequence genome;
 final String chrom_;
 final int chrom_index;
 final String[] type_nmes;
 final int num_sources;
 Sequence chr5prime; Sequence chr3prime;
 
 Stack<IdentityProfile1> idents = new Stack<IdentityProfile1>();
 
 
 public synchronized void replace(IdentityProfile1 pr1){
	 //pr1.clear();
	// System.err.print("returning");
	 idents.push(pr1);
	 
 }
 /** supp is in case this is a supplementary alignment, in which case we need to get the same profile as previously used 
  * we assume that the primary alignment is first, and supplementary alignments follow
  * note that this breaks multi-threading
  * */
 public synchronized IdentityProfile1 get(String readName, boolean supp){
	 int len = idents.size();
	// System.err.println(idents.size());
	 if(supp){
		 
		 for(int i=0; i<len; i++){
			 if(idents.get(i).readName.equals(readName)){
				 // in this case we dont clear it
				 return idents.remove(i);
			 }
		 }
		 throw new RuntimeException("we couldnt find the primary read processor - maybe this file is not sorted by read ID");
	 }
	 IdentityProfile1 idp;
	 if(len==0){
		idp =  new IdentityProfile1(this) ;
	 }else{
		 idp = idents.pop();
		 if(ViralTranscriptAnalysisCmd2.allowSuppAlignments) idp.commit(); // should be safe to commit previous
		 idp.clear(); // gets object clear for next use

	 }
	
	 if(supp && len>1){
		 // need to check how this works with multi-threading
		 throw new RuntimeException("this may not work with multi-threading");
	 }
	 idp.setName(readName);
	// System.err.println(idents.size());
	return idp;
 }
 
 public static int cleanUpMoreThan = 100000; //  distance between previous end and current start to remove cluster
	public static int clearUpThreshold = 100; //number of clusters before try to clear up 
	
	/*
 public void clearUpTo(int endThresh){
	 
	 this.all_clusters.clearUpTo(endThresh, o, genome, chrom_index);
 }*/
 public void printBreakPoints() throws IOException {
		if(bp!=null) this.bp.printBreakPoints(this.o, this.chrom_index);
		
	}
 
 public IdentityProfileHolder(ArrayList<Sequence> genomes, Sequence refSeq,Outputs o, String[] in_nmes,  boolean calcBreakpoints, int chrom_index, Annotation annot)
		 throws IOException {
	 this.genome=refSeq;
	 Arrays.fill(basesStart,0);
	 Arrays.fill(basesEnd,0);
	this.chrom_ =genome.getName();
	this.chrom_index = chrom_index;
	this.type_nmes = in_nmes;
		this.num_sources = in_nmes.length;
		this.o  = o;
	int seqlen = refSeq.length();
		all_clusters =new CigarClusters(refSeq,  num_sources, annot);
		if(calcBreakpoints){
			System.err.println("calculating break point usage");
			bp = new BreakPoints(num_sources, refSeq, genomes);
		}else{
			bp  = null;
		}
		Sequence chr = genome;
		chr5prime = TranscriptUtils.coronavirus ? chr.subSequence(0	, Math.min( primelen, chr.length())) : null;
		 chr3prime = TranscriptUtils.coronavirus ? chr.subSequence(Math.max(0, chr.length()-primelen), chr.length()) : null;
		
	}
 
 
	
	public synchronized void addBreakPoint(int source_index, int i, int br_i, int br_i1) {
	 	this.bp.addBreakPoint(source_index, i, br_i, br_i1);
		
	}
	public synchronized void addBreakPoint(int source_index, String i, int br_i, int br_i1) {
	 	this.bp.addBreakPoint(source_index, i, br_i, br_i1);
		
	}
	public void empty(){
		int len =idents.size();
		//if(len>1) {
			//throw new RuntimeException(len+" !!");
	//	}
		for(int i=0; i<len; i++){
			 IdentityProfile1 idp=idents.get(i);
			 if(ViralTranscriptAnalysisCmd2.allowSuppAlignments) idp.commit();
			idp.clear();
		}
	}
	SortedSet<String> geneNames = new TreeSet<String>();
	public void getConsensus() throws IOException {
		//waitOnThreads(100);
		this.empty();
		this.all_clusters.getConsensus( o, this.genome, this.chrom_index, this.geneNames);
		
}
	SortedMap<Integer, Integer> currentStart = new TreeMap<Integer, Integer>();
	private synchronized void removeStart(int startPos){
		int  v1 = currentStart.get(startPos);
		if(v1==0) currentStart.remove(startPos);
		else currentStart.put(startPos, v1-1);
	}
	private synchronized void addStart(int startPos){
		Integer v = currentStart.get(startPos);
		currentStart.put(startPos, v==null ? 1 : v+1);
	}
	final Integer[] basesStart = new Integer[4];final Integer[] basesEnd = new Integer[4];
	final Integer[] readCounts = new Integer[] {0,0,0};
	public void identity1(Sequence readSeq, SAMRecord sam,
			int source_index, boolean cluster_reads,  String pool, double qval, boolean supp) {
		// TODO Auto-generated method stub
		final int start = sam.getStart();
		addStart(start);
		if(ViralTranscriptAnalysisCmd2.sorted && this.all_clusters.l.size()> clearUpThreshold){ 
		 this.all_clusters.clearUpTo(currentStart.firstKey() -cleanUpMoreThan , o, genome, chrom_index);
		}
		//Boolean backwardStrand =null;
		
		Runnable run = new Runnable(){
			public void run(){
				final IdentityProfile1 profile = get(sam.getReadName(), supp);
				if(profile.coRefPositions.breaks.size()>0 && ! supp) throw new RuntimeException("this is not clear");
			//	if(supp && profile.suppl!=null && profile.suppl.size()>0) throw new RuntimeException("supps not empty");
				
				profile.identity1(genome, chr5prime,chr3prime, readSeq, sam, source_index, cluster_reads,  pool, qval, supp);
			//	profile.commit();
				replace(profile);
				removeStart(start);
			}
		};
		if(executor==null) run.run();
		else executor.execute(run);
	

	}	
	
	public static void shutDownExecutor() {
		if(executor!=null) {
			waitOnThreads(10);
			executor.shutdown();
		}
		
	}
	
 
}
