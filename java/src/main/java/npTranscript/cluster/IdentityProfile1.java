package npTranscript.cluster;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.HashSet;
import java.util.Iterator;
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
	public static boolean trainStrand = true;
	public static int min_first_last_exon_length= 20;
	/*final String[] type_nmes; 
	public  Sequence ref;
	public String chrom_;
	
	public final Outputs o;
	public CigarClusters all_clusters;
	final public BreakPoints bp;
	*/
	public  int chrom_index;
	String chrom_;
	
	
	
	public static boolean annotByBreakPosition = true;
	public static int writeCoverageDepthThresh = 100;
	public static int[] writeIsoformDepthThresh = new int[] {10};
	public static int msaDepthThresh = 10;
	public static boolean includeStartEnd = true;
	
	public static boolean attempt5rescue = false;
	public static boolean attempt3rescue = false;
	public static int extra_threshold = 500;
	public static int extra_threshold1 = 50;
	public static int extra_threshold2 = 20;
	public static double qual_thresh = 20.0D;
	public static int break_thresh = 1000;
	
	public static Sequence polyA = new Sequence(Alphabet.DNA(), "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA".toCharArray(), "polyA");
	//public  static boolean tryComplementOnExtra = false;
	public static boolean reAlignExtra = false;
	
	public int startPos, endPos;
	public int readSt, readEn; 
	
	static Integer[] seq1 = new Integer[4];
	static Integer[] seq2 = new  Integer[2];
	static Integer[] seq11 = new Integer[4];
	static Integer[] seq21 = new  Integer[2];

	
	final IdentityProfileHolder parent;
	final Outputs o; 
	final CigarClusters all_clusters;
	public IdentityProfile1(IdentityProfileHolder parent) {
		this.parent=parent;
		//this.chrom_index = parent.chrom_index;
		//this.chrom_ = parent.chrom_;
		this.num_sources = parent.num_sources;
		this.coRefPositions = new CigarCluster("", "-1",num_sources,'.');
		this.o = parent.o;
		this.all_clusters = parent.all_clusters;
		
	}
	
	public void setName(String readName, String chrom_, int chrom_index) {
		// TODO Auto-generated method stub
		this.readName = readName;
		this.chrom_ = chrom_;
		this.chrom_index = chrom_index;
	}
	/*public IdentityProfile1(Sequence refSeq,
			Outputs o, CigarClusters all_clusters, BreakPoints bp, 
			String[] in_nmes,  int startThresh, int endThresh, boolean calcBreakpoints, Sequence chrom, int chrom_index) throws IOException {
	this.ref = chrom;
	this.bp = bp;
	this.all_clusters = all_clusters;
	this.chrom_ =chrom.getName();
	this.chrom_index = chrom_index;
	this.type_nmes = in_nmes;
		this.num_sources = in_nmes.length;
		this.coRefPositions = new CigarCluster("-1",num_sources,'.');
		this.genome = refSeq;
		this.source_index = 0;		
		this.o  = o;
		refBase = 0;
		readBase = 0;
		int seqlen = refSeq.length();
		
	}*/

	

	
	
	


static char delim = '_'; // delim for second key
static char delim1 = ',';
static char delim_start ='$';

	static double break_round = 10.0;
	static String NAstring = "NA";

	public static boolean subclusterBasedOnStEnd = false;
	
	
	SortedSet<String> geneNames = new TreeSet<String>();
	
	public String[] clusterID = new String[2];
	
	//** this adds coRefPositions into the clusters
	public void commit(){
		char strand = this.coRefPositions.strand;
		int start_read = this.readSt; int end_read = this.readEn;int readLength = end_read-start_read;
		 parent.all_clusters.matchCluster(coRefPositions, this.source_index, this.num_sources,  this.chrom_,this.chrom_index,  clusterID, strand, this.readName); // this also clears current cluster
	//	int len1 = readSeq.length();
		String str = id+"\t"+clusterID[0]+"\t"+clusterID[1]+"\t"+source_index+"\t"+readLength+"\t"+start_read+"\t"+end_read+"\t"
		+type_nme+"\t"+chrom_+"\t"
		+startPos+"\t"+endPos+"\t"+strand+"\t"+coRefPositions.numIsoforms()+"\t"+(hasLeaderBreak ? 1:0)+"\t"
		+coRefPositions.getError(source_index)+"\t"+secondKeySt+"\t"+strand+"\t"+breakSt+"\t"+span_str+"\t"+geneNames.size()+"\t"+String.format("%5.3g", q_value).trim();
	//	if(trainStrand){
		//	str = str+"\t"+readSeq.subSequence(0, 10)+"\t"+readSeq.subSequence(len1-10, len1)
//				+"\t"+toString(phredQ,0,10)+"\t"+toString(phredQ,len1-10,len1)
	//			+"\t"+ViralChimericReadsAnalysisCmd.median(phredQ,0,10)+"\t"+ViralChimericReadsAnalysisCmd.median(phredQ,len1-10,10)+"\t"+annot.getStrand(coRefPositions.span.iterator());			
		//}
		parent.o.printRead(str);
		
	}
	boolean hasLeaderBreak;		String breakSt; String secondKeySt;boolean includeInConsensus = true;
	byte[] phredQ; String baseQ; String type_nme; String span_str; int span;  double q_value; String id;

	/*Note:  align5prime will be null if not coronavirus */
	public String processRefPositions( SAMRecord sam, String id, boolean cluster_reads, Sequence refSeq1, int src_index , Sequence readSeq, String baseQ, 
			byte[] phredQ,
			int start_read, int end_read, char strand, SWGAlignment align5prime, SWGAlignment align3prime,
			SWGAlignment align3primeRev,
			Integer offset_3prime, Integer polyAlen, String pool, double q_value, Annotation annot
			) throws IOException, NumberFormatException{
		//CigarHash2 breaks  = coRefPositions.breaks;
		this.baseQ  = baseQ; this.phredQ = phredQ; this.q_value = q_value;
		this.coRefPositions.strand = strand; this.id = id; 
		CigarHash2 breaks = this.coRefPositions.breaks;
		startPos = breaks.get(0);
		endPos= breaks.get(breaks.size()-1);
		//startPos = sam.getAlignmentStart();//+1; // transfer to one based
		//endPos = sam.getAlignmentEnd();
		
		readSt = start_read; readEn = end_read;
		if(ViralTranscriptAnalysisCmd2.RNA[source_index]){
			coRefPositions.forward = !sam.getReadNegativeStrandFlag();
		}else{
			coRefPositions.forward = null;
		}
		Boolean forward = coRefPositions.forward;
		//boolean hasSplice = false;
		int  readLength = readSeq.length();
		//Annotation annot =parent.all_clusters.annot;

	
		
	
	
		if(polyAlen!=null && polyAlen>0){
			Integer seqlen = refSeq1.length();
			seqlen = seqlen-(polyAlen);
			annot.adjust3UTR(seqlen);
		}
	
	//	coRefPositions.start = startPos;
	//	coRefPositions.end = endPos;
		coRefPositions.forward = forward;
		//breaks.add(coRefPositions.end);
		
		if( align5prime!=null ){
			if(align5prime.getIdentity()>0.85 * Math.max(start_read,align5prime.getLength())){
				System.err.println("rescued 5' "+readSeq.getName()+" "+align5prime.getIdentity()+" "+align5prime.getLength());
			
				int newStartPos = align5prime.getStart2() + 1; // transfer to 1-based
				int newBreakPos = newStartPos + align5prime.getSequence2().length -  align5prime.getGaps2();
				if(newBreakPos < startPos){
					readSt = align5prime.getStart1();
					startPos = newStartPos;
					//coRefPositions.start = newStartPos;
					coRefPositions.breaks.add(0, newBreakPos);
					coRefPositions.breaks.add(0, newStartPos);
				}
			}
		}
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
		if( align3prime!=null ){
			
			//System.err.println(readSeq.subSequence(end_read, readSeq.length()));
			double diff = readLength-end_read;
			double ident = align3prime.getIdentity();
			double alignLen = align3prime.getLength();
			if(ident > 0.85 *Math.max(alignLen,diff)){// && read_st< 20){
				includeInConsensus = false;
				int newBreakPos = offset_3prime+ align3prime.getStart2() + 1; // transfer to 1-based
				int newEndPos = newBreakPos + align3prime.getSequence2().length -  align3prime.getGaps2();
//				int read_end = read_st + align3prime.getSequence1().length - align3prime.getGaps1();
			//	System.err.println(newBreakPos+"-"+newEndPos);
				if(newEndPos > endPos){
					this.readEn = end_read+align3prime.getStart1()+align3prime.getLength()-align3prime.getGaps1();
					System.err.println("rescued 3' "+readSeq.getName());
					endPos = newEndPos;
				//	coRefPositions.end = newEndPos;
					coRefPositions.breaks.add( newBreakPos);
					coRefPositions.breaks.add( newEndPos);
				}
			}
		}
		
		//String roundStartP = annotByBreakPosition ? TranscriptUtils.round(startPos, CigarHash2.round)+"" : 	"";
		StringBuffer secondKey =new StringBuffer();
		
	
		 hasLeaderBreak = TranscriptUtils.coronavirus? (breaks.size()>1 &&  annot.isLeader(breaks.get(1))) : false;
		if(pool!=null) {
			String[] pool_ = pool.split(":");

			secondKey.append(pool_[0]+delim1); //+annot.nextDownstream(startPos,chrom_index, forward)+delim);
			if(pool_.length>1){
				parent.addBreakPoint(src_index, pool_[0], Integer.parseInt(pool_[1]), startPos);
			}
		}
		
		if(ViralTranscriptAnalysisCmd2.coronavirus){
			if(includeStartEnd || breaks.size()==2){
				secondKey.append(annot.nextUpstream(startPos,chrom_index, forward)+delim);
			}
			boolean firstBreak=true;
			for(int i=1; i<breaks.size()-1; i+=2){
				int gap = breaks.get(i+1)-breaks.get(i);
				String upst = annot.nextUpstream(breaks.get(i), chrom_index,forward); //5prime break
				secondKey.append(upst+delim1);
				secondKey.append(annot.nextDownstream(breaks.get(i+1), chrom_index,forward)+delim);  //3prime break
				if(gap > break_thresh){
					if(firstBreak){
						firstBreak=false;
						if(annotByBreakPosition) parent.addBreakPoint(source_index, 0, breaks.get(i), breaks.get(i+1));
					//	if(bp!=null) this.bp.addBreakPoint(source_index, 0, breaks.get(i), breaks.get(i+1));
						hasLeaderBreak = true;
					}else{
						if(annotByBreakPosition)  parent.addBreakPoint(source_index, 1, breaks.get(i), breaks.get(i+1));

					}
				}
		
			}
			if(includeStartEnd || breaks.size()==2) {
				secondKey.append(annot.nextUpstream(breaks.get(breaks.size()-1), chrom_index, forward)); //last break is upstream start pos
			}
				
		}else if(annot instanceof EmptyAnnotation){
			if(forward!=null){
				if(forward) secondKey.append(annot.nextUpstream(breaks.get(breaks.size()-1), chrom_index, forward));
				else secondKey.append(annot.nextUpstream(breaks.get(0), chrom_index, forward));
			}else{
				secondKey.append(annot.nextUpstream(breaks.get(breaks.size()-1), chrom_index, forward));
			}
		}else{
			Collection<String> genes = new HashSet<String>();
			List<String> coords = new ArrayList<String>();
			{
				for(int i=0; i<breaks.size(); i+=2){
					int l1 = breaks.get(i);
					int r1 = breaks.get(i+1);
					Integer upst = ((GFFAnnotation)annot).nextUpstream(l1, r1, chrom_index, forward);
					if(upst!=null){
						genes.add(annot.genes.get(upst));
					}else{
						String coord = ((GFFAnnotation)annot).getCoord(forward == null || forward ? breaks.get(i+1) : breaks.get(i), chrom_index,forward);
						coords.add(coord);
					}
				}
				if(genes.size()>0){
					for(Iterator<String> it = genes.iterator(); it.hasNext();){
						secondKey.append(it.next()+delim);
					}
				}else{
					secondKey.append(forward==null || forward ?  coords.get(coords.size()-1): coords.get(0));
				}
			}
		}
		
		if(Annotation.enforceStrand){
			secondKey.append(forward ? '+' : '-');
		}
		//System.err.println(secondKey);

		type_nme = ViralTranscriptAnalysisCmd2.coronavirus  ?  annot.getTypeNme( startPos, endPos, forward) : "NA"; //coRefPositions.getTypeNme(seqlen);
		geneNames.clear();
		this.coRefPositions.setStartEnd(startPos,endPos,src_index);
		 span_str = annot.getSpan(coRefPositions.breaks, forward,  coRefPositions.span, geneNames);
		
		 span = geneNames.size();
		if(!TranscriptUtils.coronavirus && span==0 && q_value < ViralTranscriptAnalysisCmd2.fail_thresh1){
			return secondKey.toString();
		}
		breakSt = coRefPositions.breaks.toString();
		//coRefPositions.breaks.adjustBreaks(annot);
		// need to group by start position if we annotating by break pos,e.g. so 5'mapping reads map together
		secondKeySt = secondKey.toString();
		coRefPositions.breaks_hash.setSecondKey(secondKeySt);//, endPos);
	//	commit();
		
		int num_exons =(int) Math.floor( (double)  coRefPositions.breaks.size()/2.0);

		boolean writeMSA = Outputs.doMSA!=null && Outputs.msa_sources !=null && includeInConsensus  && Outputs.msa_sources.containsKey(source_index);
		if(includeInConsensus && TranscriptUtils.coronavirus){
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
		}
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
		}
		
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
	
	

	


	static class CigarClusterRepo{
		final int num_sources;
		CigarClusterRepo(int num_sources){
			this.num_sources = num_sources;
		}
		Stack<CigarCluster> s = new Stack<CigarCluster>();
	
		public synchronized CigarCluster get(){
			if(s.size()==0){
				return new CigarCluster("","-1",num_sources,'.');
			}
			return s.pop();
		}
		public synchronized void replace(CigarCluster ccr){
			s.push(ccr);
		}
	}
	
	final CigarCluster coRefPositions; // this is a temporary object used to store information for each read as it comes through 
	
	
	
	
	
	public int refBase, readBase;

	//private  PrintWriter readClusters;
	private Sequence genome;
	private final int num_sources;
	public int source_index=-1; //source index
	public void updateSourceIndex(int i) {
		this.source_index = i;
			coRefPositions.readCount[source_index]=1;
		
	}
	

	 CigarHash2 suppl_read = null; //keeps track of positions on read of all alignments (primary plus suppl)
	 List<CigarHash2> suppl = null;
	//
	public void processAlignmentBlocks(SAMRecord sam, CigarCluster coRefPositions1){
		List<AlignmentBlock> li = sam.getAlignmentBlocks();
	
		for(int i=0; i<li.size(); i++){
			AlignmentBlock ab = li.get(i);
			coRefPositions1.addBlock(ab.getReferenceStart(), ab.getLength(), i==0, i==li.size()-1);
		}
		if(coRefPositions1.breaks.size() %2 !=0) {
			throw new RuntimeException("!!");
		}
	}
	
	
	
	public void processCigar(SAMRecord sam, Sequence refSeq, Sequence readSeq,CigarCluster coRefPositions1){
		IdentityProfile1 profile = this;
		int readPos = 0;// start from 0
		int refPos = sam.getAlignmentStart() - 1;// convert to 0-based index
		coRefPositions1.addStart(sam.getAlignmentStart());
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
					coRefPositions1.add(refPos+i+1, this.source_index, false);
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
						coRefPositions1.add(refPos+i+1, this.source_index, true);
						//profile.addRefPositions(refPos + i, true);
					}
					else{
					//	profile.mismatch[refPos + i]++;
						coRefPositions1.add(refPos+i+1, this.source_index, false);
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
					coRefPositions1.add(refPos+i+1, this.source_index, true);
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
					coRefPositions1.add(refPos+i+1, this.source_index, false);
					//profile.addRefPositions(refPos+i, false);
				}
				//profile.addRefPositions(refPos);
			//	profile.mismatch[refPos] += length;
				break;
			default:
				throw new IllegalStateException("Case statement didn't deal with cigar op: " + e.getOperator());
			}// case
		} // for
		coRefPositions1.addEnd(sam.getAlignmentEnd());
		if(coRefPositions1.breaks.size() %2 !=0) {
			throw new RuntimeException("!!");
		}
	}

	/**
	 * Get the identity between a read sequence from a sam and a reference sequence
	 * commit indicates immediately commit the CigarHash object
	 * @param refSeq
	 * @param sam
	 * @return
	 */
	public  void identity1(Sequence refSeq, Sequence fivePrimeRefSeq, Sequence threePrimeRefSeq,   Sequence readSeq, SAMRecord sam, 
			int source_index, boolean cluster_reads,  String poolID, double qval, boolean supplementary, Annotation annot) throws NumberFormatException{
		
		if(supplementary){
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
			 int rSt =  sam.getReadPositionAtReferencePosition(sam.getAlignmentStart());
			 int rEnd =  sam.getReadPositionAtReferencePosition(sam.getAlignmentEnd());
			 if(this.suppl_read.overlaps(rSt, rEnd-rSt+1)>=0){
				 return ;
			 }
			this.coRefPositions.breaks.clear();
		}
		CigarCluster coref = this.coRefPositions;
		//if(coref.forward==null) throw new RuntimeException("!!");
		//int seqlen = refSeq.length();
		IdentityProfile1 profile = this;
		//Outputs output = profile.o;
	
		String id = sam.getReadName();
		//String chrom = refSeq.getName();
		profile.updateSourceIndex(source_index);
		//List<AlignmentBlock> li = sam.getAlignmentBlocks();
	//	li.get(0).get
		if(CigarCluster.recordDepthByPosition || true ){
			this.processCigar(sam, refSeq, readSeq, coref);
		}else{
			this.processAlignmentBlocks(sam, coref);
		}
		try{
			int maxl = 100;
			int tol=5;
			int st_r = sam.getReadPositionAtReferencePosition(sam.getAlignmentStart());
			int end_r = sam.getReadPositionAtReferencePosition(sam.getAlignmentEnd());
			int diff_r = readSeq.length() -end_r;
			
			
			
			String baseQ = sam.getBaseQualityString();
			byte[]phredQs = sam.getBaseQualities();
			char strand = sam.getReadNegativeStrandFlag() ? '-': '+';
			boolean forward = !sam.getReadNegativeStrandFlag();
			SWGAlignment align_5prime = null;
			SWGAlignment align_3prime = null;
			SWGAlignment align_3primeRev = null;

			int offset_3prime =0;
			Integer polyAlen = TranscriptUtils.coronavirus ? TranscriptUtils.polyAlen(refSeq) : null;
			if(TranscriptUtils.coronavirus && st_r > extra_threshold1 && attempt5rescue && st_r < 1000 && sam.getAlignmentStart()> 100 ){
				//good candidate for a missed 5' alignment
				try{
				 align_5prime = SWGAlignment.align(readSeq.subSequence(0, st_r), fivePrimeRefSeq);
				}catch(Exception exc){
					System.err.println("warning could not attemp 5'rescue "+ readSeq.subSequence(0, st_r));
				}
			}
			
			Integer seqlen1 = TranscriptUtils.coronavirus ? refSeq.length()-polyAlen : null;
			if(TranscriptUtils.coronavirus && diff_r > extra_threshold2 && attempt3rescue &&  diff_r < 1000 && sam.getAlignmentEnd()< seqlen1- 100 ){
				//good candidate for a missed 3' alignment
				try{
				 align_3prime = SWGAlignment.align(readSeq.subSequence(end_r+1,readSeq.length()), threePrimeRefSeq);
				// align_3primeRev = SWGAlignment.align(TranscriptUtils.revCompl(readSeq.subSequence(end_r+1,readSeq.length())), threePrimeRefSeq);

				 offset_3prime = refSeq.length()-threePrimeRefSeq.length();
				}catch(Exception exc){
					System.err.println("warning could not attempt 3'rescue "+ readSeq.subSequence(end_r+1,readSeq.length()));
				}
			}
			Sequence leftseq = null;  Sequence rightseq = null;
		 double phredQL = st_r < 20 ? 0 : npTranscript.run.ViralChimericReadsAnalysisCmd.median(phredQs, 0,st_r);
		// double phredQ = npTranscript.run.ViralChimericReadsAnalysisCmd.median(phredQs, st_r, end_r - st_r);
		 double phredQR = diff_r <20 ? 0 : npTranscript.run.ViralChimericReadsAnalysisCmd.median(phredQs,  end_r , diff_r);
		 //System.err.println(phredQL+" "+phredQ+" "+phredQR);
			if(st_r> extra_threshold && Outputs.writePolyA && !Double.isNaN(phredQL) && phredQL >= qual_thresh){
				leftseq = readSeq.subSequence(0,st_r);
				try{
				SWGAlignment polyAlign =  SWGAlignment.align(leftseq, polyA);
				if(polyAlign.getIdentity() > 0.9 * polyAlign.getLength()  && polyAlign.getLength()>15){ // at least 15 A wih
					int st = polyAlign.getStart1();
					int end = st + polyAlign.getLength() - polyAlign.getGaps1();
					if(st> 10){
						String nme = readSeq.getName();
						profile.o.writePolyA(readSeq.subSequence(0, end),  nme+".L",
								baseQ.length()==1 ? baseQ: baseQ.substring(0,end), sam.getReadNegativeStrandFlag(),source_index);
						profile.o.writePolyA(readSeq.subSequence( end+1, readSeq.length()),  nme+".R",
								baseQ.length()==1 ? baseQ: baseQ.substring(end+1, readSeq.length()), sam.getReadNegativeStrandFlag(),source_index);
						System.err.println("internal polyA 5´ ");
						return;
					}
				}
				}catch(Exception exc){
					System.err.println("warning could not do polyA alignment");
				}

			}
			if(readSeq.length() -  end_r >extra_threshold && Outputs.writePolyA &&  !Double.isNaN(phredQR) && phredQR >= qual_thresh){
				rightseq = readSeq.subSequence(end_r+1, readSeq.length());
				try{
				SWGAlignment polyAlign =  SWGAlignment.align(rightseq, polyA);
				if(polyAlign.getIdentity() > 0.9 * polyAlign.getLength()  && polyAlign.getLength()>10 ){
					 int st = polyAlign.getStart1() + end_r;
					 int end = st + polyAlign.getLength() - polyAlign.getGaps1();
					 if(end< readSeq.length()-10){
						 String nme = readSeq.getName();
						
							profile.o.writePolyA(readSeq.subSequence(0, end), nme+".L",
									baseQ.length()==1 ? baseQ: baseQ.substring(0,end), sam.getReadNegativeStrandFlag(),source_index);
							profile.o.writePolyA(readSeq.subSequence( end, readSeq.length()),  nme+".R",
									baseQ.length()==1 ? baseQ: baseQ.substring(end, readSeq.length()), sam.getReadNegativeStrandFlag(),source_index);
							System.err.println("internal polyA 3'");
					return;
					 }
				}
				}catch(Exception exc){
					System.err.println("warning could not do polyA alignment");
				}
			}
		//	CigarHash2 breaks;
			if(supplementary){
				 this.suppl.add(coref.breaks.clone());
				 coref.breaks  =CigarHash2.merge(suppl);
				 coref.startPos = coref.breaks.get(0);
				 coref.endPos = coref.breaks.get(coref.breaks.size()-1);
			 }
			
			String secondKey= profile.processRefPositions( sam, id, cluster_reads, 
					refSeq, source_index, readSeq,baseQ, phredQs, st_r, end_r, strand, 
					align_5prime, align_3prime,align_3primeRev, offset_3prime, polyAlen, poolID, qval, annot);
			if(!ViralTranscriptAnalysisCmd2.allowSuppAlignments)	profile.commit();
			
			if( supplementary){
				 suppl_read.add(this.readSt); suppl_read.add(readEn);
				 Collections.sort(suppl_read);
			}
			//String ID = profile.clusterID
			diff_r = readSeq.length() - profile.readEn;
			seq11[0]=  profile.readSt; seq11[1] = profile.readEn; 
			seq11[2] = profile.startPos; seq11[3] = profile.endPos;
			
		//	System.err.println(refSeq.length());
			//check read not mapping to start or end of reference (?)
			if( (st_r>extra_threshold || end_r>extra_threshold)){
			//String desc = ";start="+sam.getAlignmentStart()+";end="+sam.getAlignmentEnd()+";full_length="+readSeq.length()+";strand="+strand;
				
			if(st_r >extra_threshold  && st_r > diff_r  && !Double.isNaN(phredQL) && phredQL >= qual_thresh){
				if(leftseq==null) leftseq = readSeq.subSequence(0,st_r);
				//if(leftseq==null) leftseq = readSeq.subSequence(0, st_r);
				String baseQL = baseQ.equals("*") ?  baseQ : baseQ.substring(0, st_r);

				/*if(sam.getReadNegativeStrandFlag()) {
					leftseq =TranscriptUtils.revCompl(leftseq); 
					leftseq.setDesc(leftseq.getDesc()+" rev "+chrom);		
					baseQL = (new StringBuilder(baseQ)).reverse().toString();
				}
				else{
					leftseq.setDesc(leftseq.getDesc()+" fwd "+chrom);		
				}*/
			
				leftseq.setName(readSeq.getName());
				StringBuffer desc = new StringBuffer();
				desc.append("L "+secondKey+" "+st_r+" "+getString(seq11));
			
				if(reAlignExtra){
					Sequence refSeq1 = refSeq.length() < 30000  ? refSeq : 
							refSeq.subSequence(profile.startPos - 10000, profile.endPos+1000);
				align_5prime = SWGAlignment.align(leftseq, refSeq1);
				 TranscriptUtils.getStartEnd(align_5prime, seq1, seq2, 0, 0, sam.getReadNegativeStrandFlag());
				String secondKey1 =  annot.nextUpstream(seq1[2], profile.chrom_index, forward)+";"+annot.nextDownstream(seq1[3], profile.chrom_index,forward);
				desc.append(" "+String.format("%5.3g",(double)seq2[1]/(double) seq2[0]).trim()+" "+secondKey1+" "+getString(seq1)+";"+getString(seq2));
				/*if(tryComplementOnExtra){
					align_5prime = SWGAlignment.align(TranscriptUtils.revCompl(leftseq), refSeq);
					 TranscriptUtils.getStartEnd(align_5prime, seq1, seq2, 0, 0, !sam.getReadNegativeStrandFlag());
						 secondKey1 =  profile.all_clusters.annot.nextUpstream(seq1[2], profile.chrom_index,forward)+";"+profile.all_clusters.annot.nextDownstream(seq1[3], profile.chrom_index,forward);

					 desc.append(" "+secondKey1+" "+String.format("%5.3g",(double)seq2[1]/(double) seq2[0]).trim()+" "+getString(seq1)+";"+getString(seq2));
						
				}*/
				}
				
				leftseq.setDesc(desc.toString());

				//leftseq.setDesc("len="+leftseq.length()+desc);//+";mtch5="+mtch_5+";mtch_3="+mtch_3);
				profile.o.writeLeft(leftseq,baseQL,sam.getReadNegativeStrandFlag(), source_index);// (double) Math.max(mtch_3, mtch_5)> 0.7 * (double)leftseq1.length());
				
			}
			if(diff_r >extra_threshold && diff_r> st_r && !Double.isNaN(phredQR) && phredQR >= qual_thresh){
				if(rightseq==null)			rightseq = readSeq.subSequence(end_r+1, readSeq.length());
			//	Sequence rightseq = readSeq.subSequence(end_r, readSeq.length());
			///	Sequence spanning1 = refSeq.subSequence(Math.max(0, sam.getAlignmentEnd()-10),sam.getAlignmentEnd());
			//	Sequence spanning2 = refSeq.subSequence(sam.getAlignmentEnd(),Math.min(refSeq.length(), sam.getAlignmentEnd()+10));
				
				 rightseq.setName(readSeq.getName());
					String baseQR = baseQ.equals("*") ?  baseQ : baseQ.substring(end_r, readSeq.length());
				 
				if(sam.getReadNegativeStrandFlag()) {
					rightseq = TranscriptUtils.revCompl(rightseq);
					baseQR = (new StringBuilder(baseQ)).reverse().toString();
				
				}
				StringBuffer desc = new StringBuffer();
				desc.append("R "+secondKey+" "+diff_r+" "+getString(seq11));
				if(reAlignExtra){
					
					Sequence refSeq1 = refSeq.length() < 30000  ? refSeq : 
							refSeq.subSequence(profile.startPos - 10000, profile.endPos+1000);
					try{
				 align_3prime = SWGAlignment.align(rightseq, refSeq1);
					
				 TranscriptUtils.getStartEnd(align_3prime, seq1, seq2, end_r, 0, sam.getReadNegativeStrandFlag());
					String secondKey1 =  annot.nextUpstream(seq1[2], profile.chrom_index,forward)+";"+annot.nextDownstream(seq1[3], profile.chrom_index,forward);

				 desc.append(" "+String.format("%5.3g",(double)seq2[1]/(double) seq2[0]).trim()+" "+secondKey1+" "+getString(seq1)+";"+getString(seq2));
				 /*if(tryComplementOnExtra){
						align_3prime = SWGAlignment.align(TranscriptUtils.revCompl(rightseq), refSeq);
						 TranscriptUtils.getStartEnd(align_3prime, seq1, seq2, end_r, 0, !sam.getReadNegativeStrandFlag());
							 secondKey1 =  profile.all_clusters.annot.nextUpstream(seq1[2], profile.chrom_index,forward)+";"+profile.all_clusters.annot.nextDownstream(seq1[3], profile.chrom_index,forward);

						 desc.append(" "+secondKey1+" "+String.format("%5.3g",(double)seq2[1]/(double) seq2[0]).trim()+" "+getString(seq1)+";"+getString(seq2));
							
					}*/
					}catch(Exception exc){
						System.err.println("warning 3'alignment unsuccessful");
					}
				}
				 rightseq.setDesc(desc.toString());
				 profile.o.writeLeft(rightseq,baseQR, sam.getReadNegativeStrandFlag(), source_index);///,  (double) Math.max(mtch_3, mtch_5)> 0.7 * (double)rightseq1.length());
			}
			
		}
		}catch(IOException exc){
			exc.printStackTrace();
		}

	}

	

	



	private boolean checkCompatible(CigarCluster coRefPositions2, SAMRecord sam) {
		boolean strand = sam.getReadNegativeStrandFlag();
		boolean negStrand = coRefPositions.strand=='-';
		if(negStrand != strand){
			//wrong strand
			//throw new RuntimeException("wrong strand");
			
			return false;
		}
		List<AlignmentBlock> l = sam.getAlignmentBlocks();
		for(int i=0; i<l.size(); i++){
			AlignmentBlock bl = l.get(i);
			int overlap_index = coRefPositions2.breaks.overlaps(bl.getReferenceStart(),bl.getLength());
			if(overlap_index>=0){
				return false;
			}
		}
		return true;
	}



	public void clear() {
		if(suppl!=null){
			this.suppl_read.clear();
			this.suppl.clear();
		}
		this.coRefPositions.clear();
		
	}

String readName="";

	







	

}