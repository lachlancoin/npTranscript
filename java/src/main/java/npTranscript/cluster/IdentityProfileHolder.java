package npTranscript.cluster;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.Iterator;
import java.util.List;
import java.util.SortedSet;
import java.util.Stack;
import java.util.TreeSet;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.ThreadPoolExecutor;

import htsjdk.samtools.AlignmentBlock;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.SequenceUtil;
import japsa.seq.Alphabet;
import japsa.seq.Sequence;
import japsa.tools.seq.SequenceUtils;
import jdk.internal.org.jline.terminal.spi.TerminalProvider.Stream;
import npTranscript.NW.PolyAT;
import npTranscript.run.Barcodes;
import npTranscript.run.TranscriptUtils1;
import npTranscript.run.ViralTranscriptAnalysisCmd2;

public class IdentityProfileHolder {
	
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

	public static double overlap_thresh = 0.33;
	public static double overlap_max = 20;
	
	static SAMRecord groupSam(Iterator<SAMRecord> sams, List<SAMRecord> extracted){
		SAMRecord first = sams.next();
		SAMRecord primary  = null;
		if(!first.isSecondaryOrSupplementary()) primary = first;
		
		int q1 = first.getMappingQuality();
		sams.remove();
		extracted.add(first);
		int[] overl = new int[3];
		
		for(int i=0; sams.hasNext(); i++){
			SAMRecord sam = sams.next();
			
		//	int q2 = sam.getMappingQuality();
		//	if(q2>q1 && i>1) throw new RuntimeException("!!");
			Arrays.fill(overl,0);
		//	List<AlignmentBlock> l2 = sam.getAlignmentBlocks();
			for(int j=0; j<extracted.size(); j++){
				read_overlap(extracted.get(j),sam,overl);
			}
		//	System.err.println(Arrays.asList(overl));
			if( ( overl[2] <= overlap_thresh *overl[1] && overl[2] <= overlap_thresh*overl[0] && overl[2] <= overlap_max)){
				if(!sam.isSecondaryOrSupplementary()){
					primary=sam;
				}
				extracted.add(sam);
				//System.err.println("included "+overl[0]+" ,"+overl[1]+","+overl[2]);
				sams.remove();
			}else {
			//System.err.println("excluded "+overl[0]+" ,"+overl[1]+","+overl[2]);
			}
		}
		
		Collections.sort(extracted , comp);
		
		if(false){
		for(int i=0; i<extracted.size(); i++){
			SAMRecord sr = extracted.get(i);
			System.err.println(sr.getReferenceName() + " "+sr.getReadNegativeStrandFlag());
			AlignmentBlock ab1 = sr.getAlignmentBlocks().get(0);
			AlignmentBlock ab2 = sr.getAlignmentBlocks().get(sr.getAlignmentBlocks().size()-1);
			System.err.println(ab1.getReadStart()+" - "+ab1.getReferenceStart());
			System.err.println(ab2.getReadStart()+" - "+ab2.getReferenceStart());

		}
		System.err.println("h");
		}
	//	return extracted;
		//return extracted;
		//if(flipped!=null) return flipped.intValue()==0 ? fa
		return primary;
	}
  
	static int[] getStEnd(SAMRecord s1) {
		List<AlignmentBlock> l1 = s1.getAlignmentBlocks();
		boolean neg1=s1.getReadNegativeStrandFlag();
		AlignmentBlock a1 = l1.get(0);
		AlignmentBlock a2 = l1.get(l1.size()-1);
		int st1 = a1.getReadStart();
		int end1 = a2.getReadStart()+a2.getLength();
		if(neg1) {
			st1 = s1.getReadLength()-st1;
			end1 = s1.getReadLength()-end1;
		}
		return new int[] {st1,end1};
	}

	private static void read_overlap(SAMRecord s1, SAMRecord s2, int[] overl) {
		int[] se1=getStEnd(s1);
		int[] se2=getStEnd(s2);
		int st1 = se1[0]; int end1 = se1[1];
		int st2 = se2[0]; int end2 = se2[1];
	
		int overlap = CigarHash2.overlap(st2, end2,st1 , end1);
		if(overlap>overl[2] || overl[2] ==0){
			overl[0] = end2 -st2; 
			overl[1] = end1 -st1;
			overl[2] = overlap;
		}
		
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
				
				
				String readSeq = sam_1.get(0).getReadString();//new Sequence(alph, sam_1.get(0).getReadString(), sam_1.get(0).getReadName());

				Collections.sort(sam_1 , comp_q);
				if(sam_1.get(0).isSecondaryOrSupplementary()) {
					throw new RuntimeException("first should be primary");
				}
		//		int source_index = (Integer) sam_1.get(0).getAttribute(SequenceUtils.src_tag);
				String rn= sam_1.get(0).getReadName();
				final IdentityProfile1 profile = get();
				int sze = sam_1.size();
				List<SAMRecord>sam1 = new ArrayList<SAMRecord>();
				
				int limit = 1;
				inner: for(int j=0; sam_1.size()>0 && j<limit; j++){
					sam1.clear();
					SAMRecord primary = groupSam(sam_1.iterator(), sam1);
					if(primary==null){
						System.err.println("warning not primary "+j);
						continue inner;
					}
					//boolean flip = TranscriptUtils.isFlipped(primary);
					String read_strand = (String) primary.getAttribute(PolyAT.read_strand_tag);
					
					
					byte[] bq = primary.getBaseQualities();
				//	System.err.println(sam1.size()+" "+sze);
					StringBuffer chrom = new StringBuffer();
					StringBuffer strand = new StringBuffer();//sam1.get(0).getReadNegativeStrandFlag() ? '-': '+';
					StringBuffer q_str = new StringBuffer();
					StringBuffer q_str1 = new StringBuffer();
					String chr = null;
					boolean fusion = false;
					
				//	if(flip){
				//		Collections.reverse(sam1);
				//	}
					
					Iterator<SAMRecord> it = sam1.iterator();
					while(it.hasNext()){
						boolean first = chr==null;
						SAMRecord record = it.next();
						if(read_strand==null){
							strand.append("N");
						}else{
							if(read_strand.charAt(0)=='+'){
								strand.append(record.getReadNegativeStrandFlag() ? '-': '+');
							}else{
								strand.append(record.getReadNegativeStrandFlag() ? '+': '-');
							}
						}
						q_str.append(record.getMappingQuality());
						//byte[] q = record.getBaseQualities();
					//	System.err.println(q.length+" "+record.getReadLength());
						
						int qstart = record.getReadPositionAtReferencePosition(record.getAlignmentStart());
						int qend = record.getReadPositionAtReferencePosition(record.getAlignmentEnd());
						//System.err.println("qstart-end "+qstart+" "+qend);
						q_str1.append(TranscriptUtils1.median(bq, qstart, qend-qstart));
						//if(strand1!=strand) strand = 'm';
					 	chrom.append(record.getReferenceName());
					 	if(first){
							chr = record.getReferenceName();
						}else{
							fusion = true;
						}
						if(it.hasNext()){
							chrom.append(";");
							q_str.append(";");q_str1.append(";");
						}
					}
//					System.err.println(sam1.size());
					int source_index=0;
					 profile.setName(rn, chrom.toString(),strand.toString(), q_str.toString(), q_str1.toString(),source_index, fusion);
				//	if(profile.coRefPositions.breaks.size()>0 && ! supp) throw new RuntimeException("this is not clear");
				//	if(supp && profile.suppl!=null && profile.suppl.size()>0) throw new RuntimeException("supps not empty");
					
					profile.identity1(genome,  readSeq, sam1, source_index, cluster_reads);
						
				//	profile.commit();
					
					
				//	removeStart(start);
				}
				replace(profile);
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
			Boolean forward_read=null;
			
			if(RNA  ){ // if we useFlames we ignore the strand
				forward_read = true;
			}else{
				forward_read=true;  //assume restraned
			//	boolean reverseForward = true;//true; // true gives longer UMI with polyA in it.  This could possibly be  a short cut to estimating polyA
			//	int forward_read_AT =PolyAT.assign(sam, sa, true, reverseForward, start_for_pA);
			//	int backward_read_AT =PolyAT.assign(sam, sa, false,reverseForward, start_rev_pA);// this assigns tags indicating the orientation of the original read as + or -;
			
				//if(forward_read_AT < backward_read_AT && forward_read_AT <= PolyAT.edit_thresh_AT){
				//	forward_read=true;
				//}
				//if(forward_read_AT > backward_read_AT && backward_read_AT <= PolyAT.edit_thresh_AT){
				//	forward_read=false;
				//}
			}
//			System.err.println("RNA "+RNA[source_index]);
				if(forward_read!=null){
					sam.setAttribute(PolyAT.read_strand_tag, forward_read ? "+" : "-");
				}
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
				
			//	if(forward_read==null && forward_read1!=null){
						sam.setAttribute(PolyAT.read_strand_tag, forward_read ? "+" : "-");
				//}
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
			if(Annotation.enforceStrand ){
				if(forward_read==null ){
					// cannot determine strand
					if( ViralTranscriptAnalysisCmd2.exclude_indeterminate_strand) return false;
				}
				if(forward_read!=null && RNA && !forward_read ){
					if(ViralTranscriptAnalysisCmd2.exclude_polyT_strand) {
							System.err.println ("warning RNA should always be forward strand, found polyT .. excluding ");

						return false;
					}
				}
				if(forward_read!=null  &&   !RNA){
					// flip if the original read is not forward
					int flip  = forward_read.booleanValue()!= align_reverse ? 1 : 0;
					sam.setAttribute(PolyAT.flipped_tag, flip);
				}
			}
			return true;
	}

	
 
}
