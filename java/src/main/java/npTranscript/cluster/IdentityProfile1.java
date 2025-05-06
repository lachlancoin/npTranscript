package npTranscript.cluster;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.SortedSet;
import java.util.TreeSet;

import htsjdk.samtools.AlignmentBlock;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.SAMRecord;
import japsa.seq.Sequence;
import npTranscript.NW.PolyAT;
import npTranscript.run.Barcodes;
import npTranscript.run.ViralTranscriptAnalysisCmd2;

/**
 * @author Lachlan Coin
 *
 */


public class IdentityProfile1 {
	
	
	public static List<String> types_to_include = null; //5_3 5_no3 etc
	Annotation annot = null;
	public static String annotFile  = null;
	public static final boolean checkPolyA = true;
//	public static double fusion_overlap_thresh = 0.5;
	public static String[] nmes = new String[] {"5_3", "5_no3", "no5_3", "no5_no3"};
	
	public void adjust3UTR(int seqlen2) {
		
	}
	public boolean isLeader(int prev_pos) {
		return prev_pos<100;
	}
	public int getTypeInd(int start, int end, Boolean forward, int seqlen) {
		if(start <=TranscriptUtils.startThresh) return end >= seqlen -TranscriptUtils.endThresh ? 0:1;
		else return end >= seqlen -TranscriptUtils.endThresh ? 2:3;
		//return null;
	}
	public String getTypeNme(int start, int end, Boolean forward, int seqlen) {
		return nmes[getTypeInd(start,end, forward, seqlen)];
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
	public static int writeIsoformDepthThresh = 10;
	public static int msaDepthThresh = 10;
	public static boolean includeStartEnd = true;
	
	//public static boolean attempt5rescue = false;
	//public static boolean attempt3rescue = false;
//	public static int extra_threshold = 500;
//	public static int extra_threshold1 = 50;
//	public static int extra_threshold2 = 20;
	public static double qual_thresh = 20.0D;
	public static int break_thresh = 1000;
	
	public static String[] seqs = new String[] {
			"AAAAAAAAAAAAAAAAAAAA",
//			 1st adapter 		 RTA 
			 //Top: 5' - 3'
			 "GGCTTCTTCTTGCTCTTAGGTAGTAGGTTC",
		//	 Bottom:  5' - 3'
			 "GAGGCGAGCGGTCAATTTTCCTAAGAGCAAGAAGAAGCCTTTTTTTTTT",//- 3' 

//			 2nd adapter 			 RMX 			 Top: 

			 //5' - 
			 "TGATGATGAGGGATAGACGATGGTTGTTTCTGTTGGTGCTGATATTGCTTTTTTTTTTTTTATGATGCAAGATACGCAC", //- 3' 

			 //Bottom: 
//			 5' - 3'
			 "GAGGCGAGCGGTCAATTTGCAATATCAGCACCAACAGAAACAACCATCGTCTATCCCTCATCATCAGAACCTACTA" 
			 //- 3' 

	};
	
	/*public static Sequence[] polyAs;// , polyT;
	static {
	//	char[] As = new char[20];
		//Arrays.fill(As, 'A');
		//char[] Ts = new char[20];
		//Arrays.fill(As, 'T');
		
		polyAs=  new Sequence[seqs.length];
		for(int i=0; i<polyAs.length; i++){
			polyAs[i] = new Sequence(Alphabet.DNA(), seqs[i], "adaptor."+i);
		}
		//polyT= new Sequence(Alphabet.DNA(), Ts, "polyT");
	}*/
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
		this.annot = parent.annot;
		//this.chrom_index = parent.chrom_index;
		//this.chrom_ = parent.chrom_;
		this.num_sources = parent.num_sources;
		this.coRefPositions = new CigarCluster("", "-1",num_sources,".");
		this.coRefPositions1 = new CigarCluster("", "-1",num_sources,".");

		this.o = parent.o;
		//this.all_clusters = parent.all_clusters;
		
	}
	public String strand;
	public boolean reverse;
	//public boolean flip=false;
	public void setName(String readName, String chrom_, String strand, String q_value,String q_value_base, int source_index, boolean fusion) {
		// TODO Auto-generated method stub
		this.id = readName;
	//	this.flip = flip;
		this.fusion = fusion;
	//	if(!supp) fusion = false;
		//else if(chrom_index!=this.chrom_index) fusion = true;
this.strand = strand;
this.reverse = strand.charAt(0)=='+';
this.q_value_str = q_value;
this.base_q_str = q_value_base;
this.updateSourceIndex(source_index);
		this.readName = readName;
		this.chrom_ = chrom_;// fusion ? this.chrom_+";"+chrom_ : chrom_;
		if(reverse && strand.length()>1) {
			List<String> chr1= Arrays.asList(chrom_.split(","));
			Collections.reverse(chr1);
			this.chrom_ = String.join(",", chr1.toArray(new String[0]));
			this.strand = (new StringBuilder(strand)).reverse().toString();
		}
		if(ViralTranscriptAnalysisCmd2.RNA){
			coRefPositions.forward = strand.charAt(0)=='+';
		}else{
			coRefPositions.forward = null;
		}
		//this.chrom_index = chrom_index;
		//this.seqlen = seqlen;
	}

boolean fusion = false;
	
	public static List<int[]> coords = new ArrayList<int[]>();
	public static List<String> coords_name = new ArrayList<String>();
	
	
	


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
"readID\tsource\tchrom\tstartPos\tendPos\tendStr\talign_strand\ttype\terrorRatio\tid\tbreaks_ref\tbreaks_read\tread_len\tqval\tbaseQ\tpolyA\tpolyA_fusion_site";
	// public static String polyA = "AAAAAAAAAAAAAA";
	 //public static String polyT = "TTTTTTTTTTTTTT";
	// public static int polyA_ed_dist = 1;
	// public static int checkdist = 20;
	public void commit(Sequence readSeq, SAMRecord sam) throws IOException{
		//int[] res = new int[2];
		int read_len = readSeq.length();
		if(sam.isSecondaryOrSupplementary()) throw new RuntimeException("need primary ");
		String[] res1 = new String[2];
		StringBuffer st = new StringBuffer();
		if(chrom_.contains(";") && !fusion ){
			throw new RuntimeException("!!");
		}
		//boolean hasPoly = false;
		if(fusion ){
			String seq = sam.getReadString();
			;//, polyA_ed_dist);
			
			
					int len1 = coRefPositions.start_positions.size();
			 		for(int j=1; j<len1; j++){
						int br = coRefPositions1.breaks.get(coRefPositions.start_positions.get(j));
								//st.append(PolyAT.check(seq,br));
								if(j<len1-1) st.append(",");
						}
			 		
		}else{
			st.append(".");
		}
		
		String strand = this.coRefPositions.strand;
	//	int start_read = this.readSt; int end_read = this.readEn;int readLength = end_read-start_read;
		// parent.all_clusters.matchCluster(coRefPositions, this.source_index, this.num_sources,  this.chrom_,  clusterID, strand, this.readName); // this also clears current cluster
	//	int len1 = readSeq.length();
		//Integer flipped = (Integer) sam.getAttribute(PolyAT.flipped_tag);
		Integer pt = sam.getIntegerAttribute("pt");
		String pt1 = pt==null ? "NA" :pt.toString();
		String id_ = parent.o.writeString(coRefPositions, source_index, chrom_, strand);
		parent.o.addDepthMap(id_, coRefPositions.map);
		String str = id+"\t"+source_index+"\t"
		+chrom_+"\t"
		+startPos+"\t"+endPos+"\t"+endPosStr+"\t"+strand+"\t"+this.type_nme+"\t"
		+coRefPositions.getError(source_index)+"\t"+id_+"\t"+breakSt+"\t"+breakSt1+"\t"+read_len+"\t"+q_value_str.trim()+"\t"+this.base_q_str.trim()+"\t"+
		pt1+"\t"+
		(st.toString() );
	//	if(trainStrand){f
		//	str = str+"\t"+readSeq.subSequence(0, 10)+"\t"+readSeq.subSequence(len1-10, len1)
//				+"\t"+toString(phredQ,0,10)+"\t"+toString(phredQ,len1-10,len1)
	//			+"\t"+ViralChimericReadsAnalysisCmd.median(phredQ,0,10)+"\t"+ViralChimericReadsAnalysisCmd.median(phredQ,len1-10,10)+"\t"+annot.getStrand(coRefPositions.span.iterator());			
		//}
		String str1 = PolyAT.getInfo(sam);
		String str2 = ViralTranscriptAnalysisCmd2.barcodes == null ? "": Barcodes.getInfo(sam);
		parent.o.printRead(str+"\t"+str1+"\t"+str2);
		
		
		parent.o.append(chrom_, strand, endPosStr, id_, pt);
		
		
/*		int start_target=1;
		int end_target=50;
		end_target =  29903;
		start_target = end_target -100;
		
		if(st0.length>1){
		
		st0 = st0[1].split("_");
		*/
		
		if(includeInConsensus  && Outputs.msa_sources.containsKey(source_index) && 
				(Outputs.doMSA!=null && Outputs.doMSA.contains(type_nme)  )
				&& breakSt.indexOf(";")<0
				) {
			String[] st0 = breakSt.split(";")[0].split(",");
			for(int k=0; k<st0.length; k++){
				String[] st0_1 = st0[k].split("_");
				int start_ref = Integer.parseInt(st0_1[0]);
				int end_ref = Integer.parseInt(st0_1[1]);
			for(int j=0; j<coords.size(); j++){
				int[] start_end = coords.get(j);
				int start_target=start_end[0];
				int end_target = start_end[1];
				int targ_leng = end_target-start_target;

			if(start_ref < start_target+5 && end_ref > end_target-5){
				//String[] st1 = breakSt1.split(";")[0].split(",")[1].split("_");
		//	if(true){
				int start_read = sam.getReadPositionAtReferencePosition(start_target);//Integer.parseInt(st1[0]);
				int end_read=sam.getReadPositionAtReferencePosition(end_target);
				int read_leng = end_read-start_read;
				if(start_read< end_read && read_leng < targ_leng+10 && read_leng >targ_leng-10){
						Sequence readSeq1 = readSeq.subSequence(start_read,end_read );
						String baseQ1 = baseQ.length()<=1 ? baseQ :baseQ.substring(start_read-1, end_read-1);
				//		List<Integer>read_breaks = new ArrayList<Integer>();
					//	for(int i=0; i<breaks.size(); i++){
					//		read_breaks.add(sam.getReadPositionAtReferencePosition(breaks.get(i)-1, true));
					//	}
						//need to put breaks back in 0 coords for fasta file
						int chrom_index=0;
						readSeq1.setDesc(chrom_index+" "+breakSt+" "+breakSt1+" "+(end_read-start_read)+ " "+strand+" "+source_index);
						String prefix =coords_name.get(j)+".";//TranscriptUtils.coronavirus ? "": num_exons+"_";
						String subID = "";//"_"+clusterID[1]+"_"
						String secondKeySt1 = start_target+"_"+end_target;//st0[0]+"_"+st0[1];
						parent.o.writeToCluster(prefix+secondKeySt1,subID, source_index, readSeq1, baseQ1,  readSeq.getName(), strand.charAt(0));
				}
			}
			
		   }
		}
		}
		
	}
	boolean hasLeaderBreak_;		String breakSt; String breakSt1; String secondKeySt;boolean includeInConsensus = true;
	byte[] phredQ; String baseQ; String type_nme; String span_str; int span;  String q_value_str; String id;
	String base_q_str;
	
	String endPosStr;
	
	
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
		//final Annotation annot;// = this.annot==null ? this : annot;
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
		//
	boolean 	 hasLeaderBreak = ViralTranscriptAnalysisCmd2.coronavirus? (breaks.size()>1 &&  annot.isLeader(breaks.get(0))) : false;
		
		if(ViralTranscriptAnalysisCmd2.coronavirus){
			boolean forward1 = true;//orward
			if(includeStartEnd || breaks.size()==2){
				secondKey.append(annot.nextUpstreamG(startPos, forward1)+delim);
			}
			boolean firstBreak=true;
			for(int i=1; i<breaks.size()-1; i+=2){
				int gap = breaks.get(i+1)-breaks.get(i);
				String upst = annot.nextUpstreamG(breaks.get(i),forward1); //5prime break
				secondKey.append(upst+delim1);
				String downst = annot.nextDownstreamG(breaks.get(i+1),forward1);
				secondKey.append(downst+delim);  //3prime break
				if(annotByBreakPosition && gap > break_thresh){
					if(firstBreak){
						firstBreak=false;
						if(annotByBreakPosition) parent.addBreakPoint(source_index, 0, breaks.get(i), breaks.get(i+1));
					//	if(bp!=null) this.bp.addBreakPoint(source_index, 0, breaks.get(i), breaks.get(i+1));
					//	hasLeaderBreak = true;
					}else{
						if(annotByBreakPosition)  parent.addBreakPoint(source_index, 1, breaks.get(i), breaks.get(i+1));

					}
				}
		
			}
		//	if(includeStartEnd || breaks.size()==2) {
				secondKey.append(annot.nextUpstreamG(breaks.get(breaks.size()-1),  forward1)); //last break is upstream start pos
		//	}
				
		}else if(annot instanceof EmptyAnnotation){
			//secondKey.append(chrom_+"/"+strand+"/");
	
			
			for(int jjk1 =0; jjk1<this.strand.length(); jjk1++){
				
				int	jjk = reverse ? strand.length()-(1+jjk1) : jjk1;
				
				if(jjk1>0 ) secondKey.append(";");
				int s1 = startPositions.get(jjk);
				int e1 = jjk < startPositions.size()-1 ? startPositions.get(jjk+1) : breaks.size();
				/*if(flip){
					int s2 = s1;
					int e2 = e1;
					s1 = e2;
					e1 = s2;
				}*/
				///prob could rethink this
				if(forward!=null){
					if(strand.charAt(jjk)=='+') secondKey.append(this.annot.nextDownstream(breaks.get(e1-1), forward));//TranscriptUtils.round(breaks.get(e1-1), CigarHash2.round1));
					else secondKey.append(this.annot.nextDownstream(breaks.get(s1), forward));//TranscriptUtils.round(breaks.get(s1),CigarHash2.round1));
				}else{
					secondKey.append(this.annot.nextDownstream(breaks.get(e1-1), true));//TranscriptUtils.round(breaks.get(e1-1), CigarHash2.round1));
				}
				
			}
			this.endPosStr = secondKey.toString();
			//secondKey.insert(0, chrom_+"/"+strand+"/");
		}else{
			throw new RuntimeException("!!");
		}
		
		
	//	if(Annotation.enforceStrand){
			
	//	}
		//System.err.println(secondKey);
		//if(hasLeaderBreak){
		//	System.err.println("has leader break");
		//}
	//	int seqlen = readSeq.length();
		type_nme =  getTypeNme( startPos, endPos, forward, annot.seqLen()) ; //coRefPositions.getTypeNme(seqlen);
		geneNames.clear();
		this.coRefPositions.setStartEnd(startPos,endPos,src_index);
		 span_str = "";//annot.getSpan(coRefPositions.breaks, forward,  coRefPositions.span, geneNames);
		
		 span = geneNames.size();
	/*	if(!TranscriptUtils.coronavirus && span==0 && q_value < ViralTranscriptAnalysisCmd2.fail_thresh1){
			return secondKey.toString();
		}*/
		breakSt = coRefPositions.getBreakString(strand);
		breakSt1 = coRefPositions1.getBreakString(strand);

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
		/// this is new
		
		
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



	/*private static String getString(Integer[] seq12) {
		return CigarHash2.getString(Arrays.asList(seq12));
	}*/
	
	

	


	
	
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
//res1
	
	
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
	int primary_index=-1;
		for(int j=0; j<sams.size(); j++){
			
			SAMRecord sam = sams.get(j);
			if(!sam.isSecondaryOrSupplementary()) primary_index=j;
			String id = sam.getReadName();
			Sequence refSeq =(refSeqs==null || sam.getReferenceIndex() >= refSeqs.size()) ? null : 
				refSeqs.get(sam.getReferenceIndex());
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
	//	Integer flipped = .getIntegerAttribute(PolyAT.flipped_tag);
	//	boolean flip = flipped!=null && flipped.intValue()==1;
		boolean flip1 = TranscriptUtils.isFlipped(sams.get(primary_index));
		try{
			int maxl = 100;
			int tol=5;
			String baseQ = sams.get(primary_index).getBaseQualityString();
			byte[]phredQs = sams.get(primary_index).getBaseQualities();
			String secondKey= profile.processRefPositions(  id, cluster_reads, 
					 source_index, readSeq,baseQ, phredQs);
			
			//if(!ViralTranscriptAnalysisCmd2.allowSuppAlignments){
			  if(types_to_include==null || types_to_include.contains(profile.type_nme)){

				profile.commit(readSeq, sams.get(primary_index));
			}
			
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