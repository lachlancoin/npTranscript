package npTranscript.cluster;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.Iterator;
import java.util.List;
import java.util.Map.Entry;
import java.util.SortedMap;
import java.util.SortedSet;
import java.util.Stack;
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.ThreadPoolExecutor;

import com.google.gson.Gson;

import htsjdk.samtools.AlignmentBlock;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.SequenceUtil;
import japsa.seq.Alphabet;
import japsa.seq.Sequence;
import npTranscript.NW.AlignmentParameters;
import npTranscript.NW.PolyAT;
import npTranscript.run.ViralTranscriptAnalysisCmd2;

public class IdentityProfileHolder {
	public static boolean prime3 = false;
	static boolean CHECK=true;
	public static double overlap_thresh = 0.33;
	public static double overlap_max = 20;
	public static double gap_max = 20;  // max gap on read between adjacent alignments

	static Gson gson = new Gson();
	public static void waitOnThreads(int sleep) {
		if(executor==null) return;
		if(executor instanceof ThreadPoolExecutor){
	    	while(((ThreadPoolExecutor)executor).getActiveCount()>0){
	    		try{
		    	//System.err.println("IdentityProfileHolder: awaiting completion "+((ThreadPoolExecutor)executor).getActiveCount());
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
 //final CigarClusters all_clusters;
  BreakPoints bp;
 final Outputs o;
 //Sequence genome;
 // String chrom_="";
 // int chrom_index;
 final String type_nmes;
 final int num_sources;
 //Sequence chr5prime; Sequence chr3prime;
 
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
 public synchronized IdentityProfile1 get(){
	 int len = idents.size();
	 IdentityProfile1 idp;
	 if(len==0){
		idp =  new IdentityProfile1(this) ;
	 }else{
		 idp = idents.pop();
	
		 idp.clear(); // gets object clear for next use

	 }
	return idp;
 }
 
 public static int cleanUpMoreThan = 100000; //  distance between previous end and current start to remove cluster
	public static int clearUpThreshold = 100; //number of clusters before try to clear up 
	
	/*
 public void clearUpTo(int endThresh){
	 
	 this.all_clusters.clearUpTo(endThresh, o, genome, chrom_index);
 }
 public void printBreakPoints(int chrom_index) throws IOException {
		if(bp!=null) this.bp.printBreakPoints(this.o, chrom_index);
		
	}
*/ 
 
	/*public void updateChrom(Sequence refSeq, String chrname, int chrom_index) {
		// TODO Auto-generated method stub
		if(this.chrom_.equals(chrname)) return;
		 this.genome=refSeq;
		 this.chrom_ =chrname;
		 this.chrom_index = chrom_index;
		 chr5prime = TranscriptUtils.coronavirus ? refSeq.subSequence(0	, Math.min( primelen, refSeq.length())) : null;
		chr3prime = TranscriptUtils.coronavirus ? refSeq.subSequence(Math.max(0, refSeq.length()-primelen), refSeq.length()) : null;
		if(calcBreakpoints && genomes!=null){
			System.err.println("calculating break point usage");
			bp = new BreakPoints(num_sources, refSeq, genomes);
		}else{
			bp  = null;
		}
	}*/
 
	boolean calcBreakpoints;
	ArrayList<Sequence> genomes;
	Annotation annot;
 public IdentityProfileHolder(ArrayList<Sequence> genomes, Sequence refSeq,
		 Outputs o, String in_nmes,  boolean calcBreakpoints, Annotation annot)
		 throws IOException {
	 this.genomes = genomes;
	 this.annot = annot;
	 this.calcBreakpoints = calcBreakpoints;
	 Arrays.fill(basesStart,0);
	 Arrays.fill(basesEnd,0);
	 this.type_nmes = in_nmes;
	 this.num_sources = 1;//in_nmes.length;
	 this.o  = o;
	// all_clusters =new CigarClusters( num_sources, genomes);
	 bp=null;
	 if(calcBreakpoints && genomes!=null){
			System.err.println("calculating break point usage");
			bp = new BreakPoints(num_sources, refSeq, genomes);
	 }
	}
 
 static int lengthOnRead(SAMRecord s1) {
		List<AlignmentBlock> l1 = s1.getAlignmentBlocks();
		AlignmentBlock a1 = l1.get(0);
		AlignmentBlock a2 = l1.get(l1.size()-1);
		int st1 = a1.getReadStart();
		int end1 = a2.getReadStart()+a2.getLength();
		int diff = end1 -st1;
		return diff;
	}
 static int readStart(SAMRecord s1) {
		List<AlignmentBlock> l1 = s1.getAlignmentBlocks();
		AlignmentBlock a1 = l1.get(0);
		AlignmentBlock a2 = l1.get(l1.size()-1);
		int st1 = a1.getReadStart();
		int end1 = a2.getReadStart()+a2.getLength();
		int diff = end1 -st1;
		return st1;
	}
static class SamComparator implements Comparator<SAMRecord>{
	boolean quality;
	
	public SamComparator(boolean b) {
		this.quality = b;
	}

	@Override
	public int compare(SAMRecord o1, SAMRecord o2) {
		if(quality) {
			int i1 = Boolean.compare(o1.isSecondaryOrSupplementary(), o2.isSecondaryOrSupplementary());
			if(i1!=0) return i1;
			int i2 = -1*Integer.compare(o1.getMappingQuality(), o2.getMappingQuality());
			if(i2!=0) return i2;
			
			
			return 1*Integer.compare(readStart(o1),readStart(o2));
		}
		else return  Integer.compare(o1.getAlignmentBlocks().get(0).getReadStart(),
			o2.getAlignmentBlocks().get(0).getReadStart());
	}
	
};
static Comparator comp = new SamComparator(false);
static Comparator comp_q = new SamComparator(true);

static AlignmentParameters align_p =new AlignmentParameters();
	
	
	//static Pattern patt = Pattern.compile("[ACTG]");

	private static int read_overlap(int[] se1, int[] se2) {
		int st1 = se1[0]; int end1 = se1[1];
		int st2 = se2[0]; int end2 = se2[1];
	
		int overlap = CigarHash2.overlap(st2, end2,st1 , end1);
	/*	if(overlap>overl[2] || overl[2] ==0){
			overl[0] = end2 -st2; 
			overl[1] = end1 -st1;
			overl[2] = overlap;
		}*/
		return overlap;
	}
 
	
	public synchronized void addBreakPoint(int source_index, int i, int br_i, int br_i1) {
	 	this.bp.addBreakPoint(source_index, i, br_i, br_i1);
		
	}
	public synchronized void addBreakPoint(int source_index, String i, int br_i, int br_i1) {
	 	this.bp.addBreakPoint(source_index, i, br_i, br_i1);
		
	}
	
	static boolean PRINT_MAP=false;
	
	
	
	static int[] stEnd(AlignmentBlock ab) {
		int st = ab.getReadStart();
		int end = st+ab.getLength();
		return new int[] {st,end};
	}
	static int[] stEnd1(SAMRecord sr) {
		int st = sr.getAlignmentBlocks().get(0).getReadStart();
		AlignmentBlock ab = 
				sr.getAlignmentBlocks().get(sr.getAlignmentBlocks().size()-1);
		int end = ab.getReadStart()+ab.getLength();
		return new int[] {st,end};
	}
	public void empty(){
		int len =idents.size();
		//if(len>1) {
			//throw new RuntimeException(len+" !!");
	//	}
		for(int i=0; i<len; i++){
			 IdentityProfile1 idp=idents.get(i);
			// if(ViralTranscriptAnalysisCmd2.allowSuppAlignments) idp.commit();
			idp.clear();
		}
	}
	SortedSet<String> geneNames = new TreeSet<String>();
	
		
	

	/*SortedMap<Integer, Integer> currentStart = new TreeMap<Integer, Integer>();
	private synchronized void removeStart(int startPos){
		int  v1 = currentStart.get(startPos);
		if(v1==0) currentStart.remove(startPos);
		else currentStart.put(startPos, v1-1);
	}
	private synchronized void addStart(int startPos){
		Integer v = currentStart.get(startPos);
		currentStart.put(startPos, v==null ? 1 : v+1);
	}*/
	final Integer[] basesStart = new Integer[4];final Integer[] basesEnd = new Integer[4];
	final Integer[] readCounts = new Integer[] {0,0,0};
	static Alphabet alph = Alphabet.DNA();
	public void identity1(List<SAMRecord>sam_1,
			 boolean cluster_reads,     List<Sequence> genome, final int primaryIndex
			) {
		if(primaryIndex<0) {
			System.err.println("no primary index");
			return;
		}
		if(primaryIndex>0){
	System.err.println(" problem with primary index");			
return ;		
	//throw new RuntimeException("expected at zero");
		}
		//List<SAMRecord>sam_1 = new ArrayList<SAMRecord>();
		//sam_1.addAll(Arrays.asList(sam));
		if(sam_1.size()==0){
			return ;
		}
		// TODO Auto-generated method stub
		
		//final char strand2 = strand;
	//	final int start = sam.get(0).getStart();
	//	addStart(start);
	//	if(ViralTranscriptAnalysisCmd2.sorted && (!ViralTranscriptAnalysisCmd2.sequential || num_sources==1)  && this.all_clusters.l.size()> clearUpThreshold){ 
	//	 this.all_clusters.clearUpTo(currentStart.firstKey() -cleanUpMoreThan , o);
	//	}
		//Boolean backwardStrand =null;
		
	
		
		Runnable run = new Runnable(){
			public void run(){
			
				if(sam_1.size()==0) throw new RuntimeException("!!");
				boolean include = process(sam_1.get(primaryIndex));
				if(!include){
				
					return ;
				}
				Iterator<SAMRecord> sam_2 = sam_1.iterator();
				
				MultiSAMRecord primary= new MultiSAMRecord(sam_2.next(), IdentityProfile1.break_thresh, prime3);
				//System.err.println("readnme "+primary.readname);
				while(sam_2.hasNext()) {
					primary.update(sam_2.next(), true);
				}
				final IdentityProfile1 profile = get();
				int sze = sam_1.size();
				int limit = 1;
					List<Integer> overlaps = new ArrayList<Integer>();
				boolean remove = primary.checkList();
					if(!remove) {
					int source_index=0;
					 profile.setName(primary.readname, primary.sb,source_index, primary.fusion());
					profile.identity1(genome, primary.read_str, sam_1, source_index, cluster_reads);
					}
				//	profile.commit();
					
					
				//	removeStart(start);
				replace(profile);
			}

		

			private Object getStart(SAMRecord t) {
				// TODO Auto-generated method stub
				return null;
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
	


	private static boolean process(SAMRecord sam) {
		boolean reverseForward = ViralTranscriptAnalysisCmd2.reverseForward;
	//	int source_index = (Integer) sam.getAttribute(SequenceUtils.src_tag);
		//boolean FLAMES=ViralTranscriptAnalysisCmd2.barcodes==null  ? false: ViralTranscriptAnalysisCmd2.barcodes[source_index].FLAMES_BARCODES;

		boolean RNA = ViralTranscriptAnalysisCmd2.RNA;//[source_index];
		if(sam.isSecondaryOrSupplementary()){
			throw new RuntimeException();
		}
		int[] start_end_for = new int[2];
		int[] start_end_rev = new int[2];
		int[] start_for_pA = new int[2];
		int[] start_rev_pA = new int[2];
		
		boolean align_reverse = sam.getReadNegativeStrandFlag();

			String sa = sam.getReadString();
			if(ViralTranscriptAnalysisCmd2.illumina){
				sa = sam.getReadName();
			
			}else{
				if(align_reverse){  
					// this converts read back to original orientation
					sa = SequenceUtil.reverseComplement(sa);
				}
			}
			//Boolean forward_read=true;
			
//			System.err.println("RNA "+RNA[source_index]);
				//if(forward_read!=null){
					sam.setAttribute(PolyAT.read_strand_tag, "+");// forward_read ? "+" : "-");
				//}
		/*	
			if(ViralTranscriptAnalysisCmd2.barcodes!=null){
				if(forward_read==null) return false;
				
				if(RNA) throw new RuntimeException("not currently supporting RNA barcodes (but probably could)");
				try{
					int for_min;
				if(forward_read){
					for_min = ViralTranscriptAnalysisCmd2. barcodes[source_index].assign( sam,sa, true, start_end_for);
				
				}else{
					for_min = ViralTranscriptAnalysisCmd2.barcodes[source_index].assign( sam,sa, false, start_end_rev);
				}
				if( ViralTranscriptAnalysisCmd2.exclude_reads_without_barcode && for_min> Barcodes.tolerance_barcode){
					return false;
				}
				
			
				//if(forward_read!=null && forward_read1!=null){
					if(for_min <= Barcodes.tolerance_barcode) { //forward_read.equals(forward_read1)){
						int startpA;
						int startB;
						if(forward_read){ //reverse barcode at 3'end
							 startB =start_end_for[1];
							 startpA = start_for_pA[1];
						}else{
							 startB =start_end_rev[0];
							 startpA = start_rev_pA[0];
							
						}
						int diff = startpA-startB ;
					
//						System.err.println(forward_read+" "+startB+" "+startpA+" "+diff);
						if(diff< ViralTranscriptAnalysisCmd2.max_umi && diff>ViralTranscriptAnalysisCmd2.min_umi && !FLAMES){
							String umi;
							if(forward_read){
								
								//String barcode = (String)sam.getAttribute(Barcodes.barcode_forward_tag);
								//System.err.println(barcode);
								 int read_len = sa.length();
								umi = SequenceUtil.reverseComplement(sa.substring(Math.max(0,read_len - startpA), Math.min(read_len,read_len - startB)));
								//String umi1 = SequenceUtil.reverseComplement(sa.substring(read_len - startpA, read_len - startB+1));
							//	System.err.println(umi);
								//System.err.println(umi1);
							//	System.err.println("");
							}else{
								umi =   sa.substring(startB, startpA);
							}
						//	System.err.println(forward_read+" "+startB+" "+startpA+" "+diff+" "+umi);

							sam.setAttribute(Barcodes.umi_tag, umi);
						}
						
						sam.setAttribute(Barcodes.barcode_confidence_tag, "high");
					}else{
						sam.setAttribute(Barcodes.barcode_confidence_tag, "low");
					}
			//	}else{
				//	sam.setAttribute(Barcodes.barcode_confidence_tag, "NA");
			//	}
				
				}catch(Exception exc){
					exc.printStackTrace();
				}
			}
			*/
		
			return true;
	}

	
 
}
