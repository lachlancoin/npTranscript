package npTranscript.cluster;

import java.io.IOException;
import java.util.Arrays;
import java.util.List;
import java.util.SortedSet;
import java.util.Stack;
import java.util.TreeSet;

import htsjdk.samtools.AlignmentBlock;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.SAMRecord;
import japsa.bio.np.barcode.SWGAlignment;
import japsa.seq.Alphabet;
import japsa.seq.Sequence;
import npTranscript.run.ViralTranscriptAnalysisCmd2;

/**
 * @author Lachlan Coin
 *
 */


public class IdentityProfile1 {
	public static final boolean checkPolyA = true;
//	public static double fusion_overlap_thresh = 0.5;
	public static String[] nmes = new String[] {"5_3", "5_no3", "no5_3", "no5_no3"};
	public String nextDownstream(int rightBreak, int chrom_index, Boolean forward){
		return ""+TranscriptUtils.round(rightBreak,CigarHash2.round);
	}
	public String  nextUpstream(int rightBreak, Boolean forward){
	//	secondKey.append(this.chrom_);		secondKey.append("/");
	//	secondKey.append(this.strand);secondKey.append("/");
		
		return ""+TranscriptUtils.round(rightBreak,CigarHash2.round);
	}
	public void adjust3UTR(int seqlen2) {
		
	}
	public boolean isLeader(int prev_pos) {
		return prev_pos<100;
	}
	public int getTypeInd(int start, int end, Boolean forward) {
		if(start <=TranscriptUtils.startThresh) return 0; //end >= seqlen -TranscriptUtils.endThresh ? 0:1;
		else return 2; //end >= seqlen -TranscriptUtils.endThresh ? 2:3;
		//return null;
	}
	public String getTypeNme(int start, int end, Boolean forward) {
		return nmes[getTypeInd(start,end, forward)];
	}
	
	
	public static boolean trainStrand = true;
	public static int min_first_last_exon_length= 20;
	/*final String[] type_nmes; 
	public  Sequence ref;
	public String chrom_;
	
	public final Outputs o;
	public CigarClusters all_clusters;
	final public BreakPoints bp;
	*/
	//public  int chrom_index;//, seqlen;
	String chrom_;
	
	
	
	public static boolean annotByBreakPosition = true;
	public static int writeCoverageDepthThresh = 100;
	public static int[] writeIsoformDepthThresh = new int[] {10};
	public static int msaDepthThresh = 10;
	public static boolean includeStartEnd = true;
	
	//public static boolean attempt5rescue = false;
	//public static boolean attempt3rescue = false;
//	public static int extra_threshold = 500;
//	public static int extra_threshold1 = 50;
//	public static int extra_threshold2 = 20;
	public static double qual_thresh = 20.0D;
	public static int break_thresh = 1000;
	
	public static Sequence polyA;// , polyT;
	static {
		char[] As = new char[20];
		Arrays.fill(As, 'A');
		//char[] Ts = new char[20];
		//Arrays.fill(As, 'T');
		polyA= new Sequence(Alphabet.DNA(), As, "polyA");
		//polyT= new Sequence(Alphabet.DNA(), Ts, "polyT");
	}
	//public  static boolean tryComplementOnExtra = false;
	//public static boolean reAlignExtra = false;
	
	public int startPos, endPos;
	//public int readSt, readEn; 
	
	static Integer[] seq1 = new Integer[4];
	static Integer[] seq2 = new  Integer[2];
	//static Integer[] seq11 = new Integer[4];
	static Integer[] seq21 = new  Integer[2];

	
	final IdentityProfileHolder parent;
	final Outputs o; 
	//final CigarClusters all_clusters;
	public IdentityProfile1(IdentityProfileHolder parent) {
		this.parent=parent;
		//this.chrom_index = parent.chrom_index;
		//this.chrom_ = parent.chrom_;
		this.num_sources = parent.num_sources;
		this.coRefPositions = new CigarCluster("", "-1",num_sources,".");
		this.coRefPositions1 = new CigarCluster("", "-1",num_sources,".");

		this.o = parent.o;
		//this.all_clusters = parent.all_clusters;
		
	}
	public String strand;
	public void setName(String readName, String chrom_, String strand, String q_value,String q_value_base, int source_index, boolean fusion) {
		// TODO Auto-generated method stub
		this.id = readName;
		this.fusion = fusion;
	//	if(!supp) fusion = false;
		//else if(chrom_index!=this.chrom_index) fusion = true;
this.strand = strand;
this.q_value_str = q_value;
this.base_q_str = q_value_base;
this.updateSourceIndex(source_index);
		this.readName = readName;
		this.chrom_ = chrom_;// fusion ? this.chrom_+";"+chrom_ : chrom_;
		
		if(ViralTranscriptAnalysisCmd2.RNA[source_index]){
			coRefPositions.forward = strand.charAt(0)=='+';
		}else{
			coRefPositions.forward = null;
		}
		//this.chrom_index = chrom_index;
		//this.seqlen = seqlen;
	}

boolean fusion = false;
	
	
	


static char delim = '_'; // delim for second key
static char delim1 = ',';
static char delim_start ='$';

	static double break_round = 10.0;
	static String NAstring = "NA";

	public static boolean subclusterBasedOnStEnd = false;
	
	
	SortedSet<String> geneNames = new TreeSet<String>();
	
	//public String[] clusterID = new String[2];
	
	//** this adds coRefPositions into the clusters
	 public static String header = 
"readID\tsource\tchrom\tstartPos\tendPos\tstrand\terrorRatio\tid\tbreaks\tbreaks1\tqval\tbaseQ\tpolyA";

	public void commit(Sequence readSeq){
		int[] res = new int[2];
		StringBuffer st = new StringBuffer("_");
		
		if(fusion  && checkPolyA){
		//	if (this.coRefPositions.forward) {
		 		for(int j=1; j<coRefPositions.start_positions.size(); j++){
					int br = coRefPositions1.breaks.get(coRefPositions.start_positions.get(j));
			//	this.coRefPositions.start_positions.get(0);
				  boolean hasPoly = getPolyA(readSeq, polyA, Math.max(0, br-30), Math.min(readSeq.length(), br+30), res);
				  if(hasPoly){
					st.append(readSeq.subSequence(res[0], res[1]));
					st.append("_");
				  }
				}
			//}
		}
		
		
		String strand = this.coRefPositions.strand;
	//	int start_read = this.readSt; int end_read = this.readEn;int readLength = end_read-start_read;
		// parent.all_clusters.matchCluster(coRefPositions, this.source_index, this.num_sources,  this.chrom_,  clusterID, strand, this.readName); // this also clears current cluster
	//	int len1 = readSeq.length();
		
		String id_ = parent.o.writeString(coRefPositions, source_index, chrom_);
		parent.o.addDepthMap(id_, coRefPositions.map);
		String str = id+"\t"+source_index+"\t"
		+chrom_+"\t"
		+startPos+"\t"+endPos+"\t"+strand+"\t"
		+coRefPositions.getError(source_index)+"\t"+id_+"\t"+breakSt+"\t"+breakSt1+"\t"+q_value_str.trim()+"\t"+this.base_q_str.trim()+"\t"+st.toString();
	//	if(trainStrand){f
		//	str = str+"\t"+readSeq.subSequence(0, 10)+"\t"+readSeq.subSequence(len1-10, len1)
//				+"\t"+toString(phredQ,0,10)+"\t"+toString(phredQ,len1-10,len1)
	//			+"\t"+ViralChimericReadsAnalysisCmd.median(phredQ,0,10)+"\t"+ViralChimericReadsAnalysisCmd.median(phredQ,len1-10,10)+"\t"+annot.getStrand(coRefPositions.span.iterator());			
		//}
		parent.o.printRead(str);
		
	}
	boolean hasLeaderBreak_;		String breakSt; String breakSt1; String secondKeySt;boolean includeInConsensus = true;
	byte[] phredQ; String baseQ; String type_nme; String span_str; int span;  String q_value_str; String id;
	String base_q_str;
	
	
	/*Note:  align5prime will be null if not coronavirus */
	public String processRefPositions(  String id, boolean cluster_reads, 
			 int src_index , Sequence readSeq, String baseQ, 
			byte[] phredQ
		//	int start_read, int end_read, char strand, 
		//	Integer offset_3prime, 
			//Integer polyAlen, 
			// double q_value//, Annotation annot
			) throws IOException, NumberFormatException{
		//CigarHash2 breaks  = coRefPositions.breaks;
		IdentityProfile1  annot = this;
	//	int chrom_index = sam.getReferenceIndex();
		this.baseQ  = baseQ; this.phredQ = phredQ;// this.q_value = q_value;
		this.coRefPositions.strand = strand; 
		this.id = id; 
		CigarHash2 breaks = this.coRefPositions.breaks;
		List<Integer> startPositions = this.coRefPositions.start_positions;
		startPos = breaks.get(0);
		endPos= breaks.get(breaks.size()-1);
		//startPos = sam.getAlignmentStart();//+1; // transfer to one based
		//endPos = sam.getAlignmentEnd();
		
	//	readSt = start_read; readEn = end_read;
		
		Boolean forward = coRefPositions.forward;
		//boolean hasSplice = false;
		//int  readLength = readSeq.length();
		//Annotation annot =parent.all_clusters.annot;

	
	//	coRefPositions.start = startPos;
	//	coRefPositions.end = endPos;
		coRefPositions.forward = forward;
		//breaks.add(coRefPositions.end);
		
	
		if(min_first_last_exon_length>0 && breaks.size()>2){ 
			// this removes very small first or last exons
			// only apply to first and last exon of multi-exon gene
			int sze = breaks.size();
			int diff0 = breaks.get(1) - startPos;
			int diff1 = endPos - breaks.get(breaks.size()-2);
			if(sze==4  && diff0 < min_first_last_exon_length && diff1 < min_first_last_exon_length){ // make sure at least one exon left
				if(diff1<diff0){
					breaks.remove(sze-1); breaks.remove(sze-2);
				}else{
					breaks.remove(1); breaks.remove(0);
				}
			}
			else{
				if(diff1<min_first_last_exon_length){
					breaks.remove(sze-1); breaks.remove(sze-2);
				}
				if(diff0<min_first_last_exon_length){
					breaks.remove(1); breaks.remove(0);
				}
			}
			endPos = breaks.get(breaks.size()-1);
			startPos = breaks.get(0);
		}
		
		//String roundStartP = annotByBreakPosition ? TranscriptUtils.round(startPos, CigarHash2.round)+"" : 	"";
		StringBuffer secondKey =new StringBuffer();
		
	//	 hasLeaderBreak = TranscriptUtils.coronavirus? (breaks.size()>1 &&  annot.isLeader(breaks.get(1))) : false;
		
		
		if(true) {//annot instanceof EmptyAnnotation){
			secondKey.append(chrom_+"/"+strand+"/");
			for(int jjk =0; jjk<this.strand.length(); jjk++){
				if(jjk>0) secondKey.append(";");
				int s1 = startPositions.get(jjk);
				int e1 = jjk < startPositions.size()-1 ? startPositions.get(jjk+1) : breaks.size();
				if(forward!=null){
					if(strand.charAt(jjk)=='+') secondKey.append(TranscriptUtils.round(breaks.get(e1-1), CigarHash2.round1));
					else secondKey.append(TranscriptUtils.round(breaks.get(s1),CigarHash2.round1));
				}else{
					secondKey.append(TranscriptUtils.round(breaks.get(e1-1), CigarHash2.round1));
				}
				
			}
		}
		
		
	//	if(Annotation.enforceStrand){
			
	//	}
		//System.err.println(secondKey);

		type_nme =  getTypeNme( startPos, endPos, forward) ; //coRefPositions.getTypeNme(seqlen);
		geneNames.clear();
		this.coRefPositions.setStartEnd(startPos,endPos,src_index);
		 span_str = "";//annot.getSpan(coRefPositions.breaks, forward,  coRefPositions.span, geneNames);
		
		 span = geneNames.size();
	/*	if(!TranscriptUtils.coronavirus && span==0 && q_value < ViralTranscriptAnalysisCmd2.fail_thresh1){
			return secondKey.toString();
		}*/
		breakSt = coRefPositions.getBreakString();
		breakSt1 = coRefPositions1.getBreakString();

		//coRefPositions.breaks.adjustBreaks(annot);
		// need to group by start position if we annotating by break pos,e.g. so 5'mapping reads map together
		secondKeySt = secondKey.toString();
		coRefPositions.breaks_hash.setSecondKey(secondKeySt);//, endPos);
	//	commit();
		
		int num_exons =(int) Math.floor( (double)  coRefPositions.breaks.size()/2.0);

		boolean writeMSA = Outputs.doMSA!=null && Outputs.msa_sources !=null && includeInConsensus  && Outputs.msa_sources.containsKey(source_index);
		/*if(includeInConsensus && TranscriptUtils.coronavirus){
			//int st1 = startPos; //position>0 ? position : startPos; // start after break
			int st1 = breaks.get(breaks.size()-2); // start of last segment
			inner: for(int i=annot.start.size()-1; i>=0; i--){
				if( annot.start.get(i)<st1-10) break inner;
				if(endPos > annot.start.get(i)+100){
					annot.addCount(i,src_index, startPos<=TranscriptUtils.startThresh);
				}else{
					//System.err.println(annot.end.get(i)+ " "+endPos);
				}
				if(writeMSA && endPos > annot.end.get(i)){
					boolean include = false;
					inner1: for(int k=0; k<breaks.size(); k+=2){
						int st = breaks.get(k);
						int end = breaks.get(k+1);
						if(st-5 < annot.start.get(i) && end+5 > annot.end.get(i)){
							include=true;
							break inner1;
						}
					}
					if(include){ // this means read covers ORF without break
					//	annot.addCount(i,src_index, startPos<=TranscriptUtils.startThresh);
						int start_ref = annot.start.get(i)-1; // need to put back in zero coords
						int end_ref = annot.end.get(i)-1;// need to put back in zero coords
						int start_read1 =  sam.getReadPositionAtReferencePosition(start_ref, true);
						int end_read1 =  sam.getReadPositionAtReferencePosition(end_ref, true);
						double diff = end_ref - start_ref;
						double diff1 = end_read1 - start_read1;
						if(diff1>0 && diff1> 0.8 *diff){
							Sequence readSeq1 = readSeq.subSequence(start_read1, end_read1);
							String baseQ1 = baseQ.length()<=1 ? baseQ : baseQ.substring(start_read1, end_read1);
							readSeq1.setDesc(chrom_index+" "+start_ref+","+end_ref+" "+start_read1+","+end_read1+" "+(end_read1-start_read1)+" "+strand+" "+source_index);
							parent.o.writeToCluster("annot_"+annot.genes.get(i),null, source_index, readSeq1, baseQ1, readSeq.getName(), strand);
						}
					}
				}
			}
		}*/
		/*
		if(includeInConsensus  && Outputs.msa_sources.containsKey(source_index) && 
				(Outputs.doMSA!=null && Outputs.doMSA.contains(type_nme)  || 
						Outputs.doMSA!=null && Outputs.doMSA.contains(span) && Outputs.numExonsMSA.contains(num_exons) )
				) {
			Sequence readSeq1 = readSeq.subSequence(start_read, end_read);
			String baseQ1 = baseQ.length()<=1 ? baseQ :baseQ.substring(start_read, end_read);
			List<Integer>read_breaks = new ArrayList<Integer>();
			for(int i=0; i<breaks.size(); i++){
				read_breaks.add(sam.getReadPositionAtReferencePosition(breaks.get(i)-1, true));
			}
			//need to put breaks back in 0 coords for fasta file
			readSeq1.setDesc(chrom_index+" "+CigarHash2.getString(breaks,-1)+" "+CigarHash2.getString(read_breaks)+" "+(end_read-start_read)+ " "+strand+" "+source_index);
			String prefix = TranscriptUtils.coronavirus ? "": num_exons+"_";
			parent.o.writeToCluster(prefix+secondKeySt,"_"+clusterID[1]+"_", source_index, readSeq1, baseQ1,  readSeq.getName(), strand);
		}*/
		
		return secondKeySt+" "+span_str;
	}
	
	

	private String toString(byte[] phredQ, int i, int len1) {
		StringBuffer sb = new StringBuffer();
		for(int j=i; j<len1; j++){
			if(j>i)sb.append(",");
			sb.append((int) phredQ[j]);;
		}
		return sb.toString();
	}



	private static String getString(Integer[] seq12) {
		return CigarHash2.getString(Arrays.asList(seq12));
	}
	
	

	


	
	
	final CigarCluster coRefPositions; // this is a temporary object used to store information for each read as it comes through 
	final CigarCluster coRefPositions1;
	
	
	
	
	public int refBase, readBase;

	//private  PrintWriter readClusters;
	private Sequence genome;
	private final int num_sources;
	public int source_index=-1; //source index
	public void updateSourceIndex(int i) {
		this.source_index = i;
			coRefPositions.readCount[source_index]=1;
		
	}
	

	// CigarHash2 suppl_read = null; //keeps track of positions on read of all alignments (primary plus suppl)
//	 List<CigarHash2> suppl = null;
	//
	public void processAlignmentBlocks(SAMRecord sam, CigarCluster coRefPositions1,CigarCluster coRefPositionsRead){
		List<AlignmentBlock> li = sam.getAlignmentBlocks();
		for(int i=0; i<li.size(); i++){
			AlignmentBlock ab = li.get(i);
			int gap = coRefPositions1.addBlock(ab.getReferenceStart(), ab.getLength(), i==0, i==li.size()-1, sam.getReferenceIndex(), sam.getReadNegativeStrandFlag());
			coRefPositionsRead.addBlock(ab.getReadStart(), ab.getLength(), i==0, i==li.size()-1, 0, gap, false);

		}
		if(coRefPositions1.breaks.size() %2 !=0) {
			throw new RuntimeException("!!");
		}
	}
	
	
	
	public void processCigar(SAMRecord sam, Sequence refSeq, Sequence readSeq,CigarCluster coRefPositions1){
		IdentityProfile1 profile = this;
		int ref_ind = sam.getReferenceIndex();
		int readPos = 0;// start from 0
		int refPos = sam.getAlignmentStart() - 1;// convert to 0-based index
		coRefPositions1.addStart(sam.getAlignmentStart(), ref_ind);
		boolean nullchr = refSeq==null || refSeq.length()==0;
		for (final CigarElement e : sam.getCigar().getCigarElements()) {
			final int length = e.getLength();

			switch (e.getOperator()) {
			case H:
				// nothing todo
			//	profile.readClipped += length;
				break; // ignore hard clips
			case P:
			//	profile.readClipped += length;
				// pad is a kind of hard clipped ??
				break; // ignore pads
			case S:
				// advance on the reference
				//profile.readClipped += length;
				readPos += length;
				break; // soft clip read bases
			case N:
				// System.err.println(length);
				refPos += length;
				/*for (int i = 0; i < length && refPos + i < refSeq.length(); i++) {
					profile.refClipped[refPos + i] += 1;
				}*/
				// profile.refClipped += length;
				break; // reference skip

			case D:// deletion
				refPos += length;
				profile.refBase += length;
				for (int i = 0; i < length && (nullchr || refPos + i < refSeq.length()); i++) {
				//	profile.baseDel[refPos + i] += 1;
					coRefPositions1.add(refPos+i+1, this.source_index, false, ref_ind);
					//profile.addRefPositions(refPos + i, false);
				}
				//profile.numDel++;
				break;

			case I:
				readPos += length;
				profile.readBase += length;

			//	profile.baseIns[refPos] += length;
			//	profile.numIns++;
				break;
			case M:
				for (int i = 0; i < length && (nullchr || refPos + i < refSeq.length()); i++) {
					
					if (nullchr || refSeq.getBase(refPos + i) == readSeq.getBase(readPos + i)){
				//		profile.match[refPos + i]++;
						coRefPositions1.add(refPos+i+1, this.source_index, true, ref_ind);
						//profile.addRefPositions(refPos + i, true);
					}
					else{
					//	profile.mismatch[refPos + i]++;
						coRefPositions1.add(refPos+i+1, this.source_index, false, ref_ind);
						//profile.addRefPositions(refPos + i, false);
					}
				}
				profile.readBase += length;
				profile.refBase += length;

				readPos += length;
				refPos += length;
				break;

			case EQ:
				readPos += length;
				refPos += length;
				for (int i = 0; i < length && (nullchr || (refPos + i) < refSeq.length()); i++) {
					coRefPositions1.add(refPos+i+1, this.source_index, true, ref_ind);
//					profile.addRefPositions(refPos+i, true);
				}
				profile.readBase += length;
				profile.refBase += length;
				//profile.match[refPos] += length;
				break;

			case X:
				readPos += length;
				refPos += length;

				profile.readBase += length;
				profile.refBase += length;
				for (int i = 0; i < length && (nullchr || (refPos + i) < refSeq.length()); i++) {
					coRefPositions1.add(refPos+i+1, this.source_index, false, ref_ind);
					//profile.addRefPositions(refPos+i, false);
				}
				//profile.addRefPositions(refPos);
			//	profile.mismatch[refPos] += length;
				break;
			default:
				throw new IllegalStateException("Case statement didn't deal with cigar op: " + e.getOperator());
			}// case
		} // for
		coRefPositions1.addEnd(sam.getAlignmentEnd(), ref_ind);
		if(coRefPositions1.breaks.size() %2 !=0) {
			throw new RuntimeException("!!");
		}
	}

	static boolean getPolyA(Sequence readSeq, Sequence polyA, int st, int end, int[] res){
		boolean haspolyA = false;
		Sequence r1 = readSeq.subSequence(st, end);
		if(true){
		SWGAlignment polyAlign =  SWGAlignment.align(r1, polyA);
		if(polyAlign.getIdentity() > 0.85 * polyAlign.getLength()  && polyAlign.getLength()>10){ // at least 15 A wih
			res[0] = polyAlign.getStart1()+st;
			res[1]  = res[0] + polyAlign.getLength() - polyAlign.getGaps1();
			haspolyA = true;
		}
		}
		if(haspolyA)		System.err.println("found polyA "+r1);
		return haspolyA;
	}
	
	/**
	 * Get the identity between a read sequence from a sam and a reference sequence
	 * commit indicates immediately commit the CigarHash object
	 * @param refSeq
	 * @param sam
	 * @return
	 */
	public  void identity1(List<Sequence> refSeqs,   Sequence readSeq, List<SAMRecord> sams, 
			int source_index, boolean cluster_reads) throws NumberFormatException{
		
		
		
		
		
		IdentityProfile1 annot = this;
		/*if(supplementary){
			if(!checkCompatible(this.coRefPositions, sam)) return;
			 if(suppl==null){
				 suppl = new ArrayList<CigarHash2>();
				 this.suppl_read = new CigarHash2(false);
			 }
			 if(suppl.size()==0){
				 suppl.add(this.coRefPositions.breaks.clone()); // need to add the primary alignment
				 suppl_read.add(this.readSt); suppl_read.add(readEn);
				 Collections.sort(suppl_read);
			 }
			// List<AlignmentBlock>l = sam.getAlignmentBlocks();
			 int rSt =  sam.getReadPositionAtReferencePosition(sam.getAlignmentStart());
			 int rEnd =  sam.getReadPositionAtReferencePosition(sam.getAlignmentEnd());
			 int i = suppl_read.overlaps(rSt, rEnd-rSt+1);
			 if(i>=0){
				 int st1 = suppl_read.get(i);
				 int end1 = suppl_read.get(i+1);
				 double overlap = suppl_read.overlap(rSt, rEnd, st1, end1);
				 double overlap_thresh = fusion ? fusion_overlap_thresh : 0; 
				 System.err.println("overlap is "+overlap+" vs "+(rEnd-rSt)+" "+(end1-st1));
				 if(overlap>=overlap_thresh * (rEnd - rSt) || overlap >= overlap_thresh * (end1-st1)){
					 return ;
				 }else{
					 System.err.println("overlap is ok");
				 }
			 }
			this.coRefPositions.breaks.clear();
		}*/
		CigarCluster coref = this.coRefPositions;
		//if(coref.forward==null) throw new RuntimeException("!!");
		//int seqlen = refSeq.length();
		IdentityProfile1 profile = this;
		//Outputs output = profile.o;
	
		for(int j=0; j<sams.size(); j++){
			
		SAMRecord sam = sams.get(j);
		String id = sam.getReadName();
		Sequence refSeq = refSeqs.get(sam.getReferenceIndex());
		//String chrom = refSeq.getName();
	//	profile.updateSourceIndex(source_index);
		//List<AlignmentBlock> li = sam.getAlignmentBlocks();
	//	li.get(0).get
		if(CigarCluster.recordDepthByPosition  ){
			this.processCigar(sam, refSeq, readSeq, coref);
		}else{
			this.processAlignmentBlocks(sam, coref, coRefPositions1);
		}
		}
		try{
			int maxl = 100;
			int tol=5;
		//	int st_r = sam.getReadPositionAtReferencePosition(sam.getAlignmentStart());
		//	int end_r = sam.getReadPositionAtReferencePosition(sam.getAlignmentEnd());
		//	int diff_r = readSeq.length() -end_r;
			String baseQ = sams.get(0).getBaseQualityString();
			byte[]phredQs = sams.get(0).getBaseQualities();
		//	char strand = sam.getReadNegativeStrandFlag() ? '-': '+';
		//	boolean forward = !sam.getReadNegativeStrandFlag();
			//int offset_3prime =0;
		
			
			String secondKey= profile.processRefPositions(  id, cluster_reads, 
					 source_index, readSeq,baseQ, phredQs);
			//if(!ViralTranscriptAnalysisCmd2.allowSuppAlignments){
				profile.commit(readSeq);
		//	}
			
		/*	if( supplementary){
				 suppl_read.add(this.readSt); suppl_read.add(readEn);
				 Collections.sort(suppl_read);
			}*/
			//String ID = profile.clusterID
		//	diff_r = readSeq.length() - profile.readEn;
			//seq11[0]=  profile.readSt; seq11[1] = profile.readEn; 
			//seq11[2] = profile.startPos; seq11[3] = profile.endPos;
			
		}catch(IOException exc){
			exc.printStackTrace();
		}

	}

	

	private static boolean checkCompatible(CigarCluster coRefPositions2, SAMRecord sam) {
		boolean strand = sam.getReadNegativeStrandFlag();
	//	boolean negStrand = coRefPositions.strand=='-';
		
		//if(!this.fusion){
		List<AlignmentBlock> l = sam.getAlignmentBlocks();
		for(int i=0; i<l.size(); i++){
			AlignmentBlock bl = l.get(i);
			int overlap_index = coRefPositions2.breaks.overlaps(bl.getReferenceStart(),bl.getLength());
			if(overlap_index>=0){
				return false;
			}
		}
		//}
		return true;
	}



	public void clear() {
	/*	if(suppl!=null){
			this.suppl_read.clear();
			this.suppl.clear();
		}*/
		this.coRefPositions.clear();
		this.coRefPositions1.clear();
		
	}

String readName="";

	







	

}