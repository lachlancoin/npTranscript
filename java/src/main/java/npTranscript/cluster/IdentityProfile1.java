package npTranscript.cluster;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.SortedSet;
import java.util.Stack;
import java.util.TreeSet;
import java.util.concurrent.ExecutorService;

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
	
//	public static  ExecutorService executor ;
	/*final String[] type_nmes; 
	public  Sequence ref;
	public String chrom_;
	
	public final Outputs o;
	public CigarClusters all_clusters;
	final public BreakPoints bp;
	*/
	final public  int chrom_index;
	final String chrom_;
	
	public static boolean annotByBreakPosition = true;
	public static int writeCoverageDepthThresh = 100;
	public static int[] writeIsoformDepthThresh = new int[] {10};
	public static int msaDepthThresh = 10;
	public static boolean includeStart = true;
	
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
		this.chrom_index = parent.chrom_index;
		this.chrom_ = parent.chrom_;
		this.num_sources = parent.num_sources;
		this.coRefPositions = new CigarCluster("-1",num_sources,'.');
		this.o = parent.o;
		this.all_clusters = parent.all_clusters;
		
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

	static double break_round = 10.0;
	static String NAstring = "NA";

	public static boolean subclusterBasedOnStEnd = false;
	
	
	SortedSet<String> geneNames = new TreeSet<String>();
	
	public String[] clusterID = new String[2];
	public String processRefPositions(SAMRecord sam, String id, boolean cluster_reads, Sequence refSeq, int src_index , Sequence readSeq, String baseQ, 
			int start_read, int end_read, char strand, SWGAlignment align5prime, SWGAlignment align3prime,
			SWGAlignment align3primeRev,
			int offset_3prime, int polyAlen, String pool, double q_value
			) throws IOException, NumberFormatException{
		
		startPos = sam.getAlignmentStart()+1; // transfer to one based
		endPos = sam.getAlignmentEnd()+1;
		
		readSt = start_read; readEn = end_read;
		boolean forward = !sam.getReadNegativeStrandFlag();
		//boolean hasSplice = false;
		int  readLength = readSeq.length();
		Annotation annot =parent.all_clusters.annot;

		CigarHash2 breaks  = coRefPositions.breaks;
		int seqlen = refSeq.length();
	
		if(polyAlen>0){
			seqlen = seqlen-(polyAlen);
			annot.adjust3UTR(seqlen);
		}
	
		coRefPositions.start = startPos;
		coRefPositions.end = endPos;
		coRefPositions.forward = forward;
		breaks.add(coRefPositions.end);
		boolean includeInConsensus = true;
		if( align5prime!=null ){
			if(align5prime.getIdentity()>0.85 * Math.max(start_read,align5prime.getLength())){
				System.err.println("rescued 5' "+readSeq.getName()+" "+align5prime.getIdentity()+" "+align5prime.getLength());
			
				int newStartPos = align5prime.getStart2() + 1; // transfer to 1-based
				int newBreakPos = newStartPos + align5prime.getSequence2().length -  align5prime.getGaps2();
				if(newBreakPos < startPos){
					readSt = align5prime.getStart1();
					startPos = newStartPos;
					coRefPositions.start = newStartPos;
					coRefPositions.breaks.add(0, newBreakPos);
					coRefPositions.breaks.add(0, newStartPos);
				}
			}
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
					coRefPositions.end = newEndPos;
					coRefPositions.breaks.add( newBreakPos);
					coRefPositions.breaks.add( newEndPos);
				}
			}
		}
		
		//String roundStartP = annotByBreakPosition ? TranscriptUtils.round(startPos, CigarHash2.round)+"" : 	"";
		StringBuffer secondKey =new StringBuffer();
		
		secondKey.append(pool);
		//String upstream, upstream2, downstream, downstream2;
		int maxg = 0;
		int maxg_ind = -1;
		int maxg2=0;
		int maxg_ind2 =-1; 
		boolean hasLeaderBreak = TranscriptUtils.coronavirus? (breaks.size()>1 &&  annot.isLeader(breaks.get(1))) : false;
		if(includeStart){
			secondKey.append(annot.nextUpstream(startPos,chrom_index, forward)+delim);
		}
		if(annotByBreakPosition){
			boolean firstBreak=true;
			for(int i=1; i<breaks.size()-1; i+=2){
				int gap = breaks.get(i+1)-breaks.get(i);
				String upst = annot.nextUpstream(breaks.get(i), chrom_index,forward);
				secondKey.append(upst+delim1);
				secondKey.append(annot.nextDownstream(breaks.get(i+1), chrom_index,forward)+delim);
				if(gap > break_thresh){
					if(firstBreak){
						firstBreak=false;
						parent.addBreakPoint(source_index, 0, breaks.get(i), breaks.get(i+1));
					//	if(bp!=null) this.bp.addBreakPoint(source_index, 0, breaks.get(i), breaks.get(i+1));
						hasLeaderBreak = true;
					}else{
						 parent.addBreakPoint(source_index, 1, breaks.get(i), breaks.get(i+1));

					}
				}
				
				if(gap>maxg && annot.isLeader(breaks.get(i))){
					maxg2 = maxg;
					maxg_ind2 = maxg_ind;
					maxg = gap;
					maxg_ind = i;
				}else if(gap>maxg2){
					maxg2 = gap;
					maxg_ind2 = i;
				}
			}
		}
		secondKey.append(annot.nextDownstream(breaks.get(breaks.size()-1), chrom_index, forward));
		
		if(Annotation.enforceStrand){
			secondKey.append(forward ? '+' : '-');
		}
		//System.err.println(secondKey);

		String type_nme = annot.getTypeNme( startPos, endPos, forward); //coRefPositions.getTypeNme(seqlen);
		geneNames.clear();
		
			
			this.coRefPositions.setStartEnd(startPos,endPos,src_index);

		String span_str = annot.getSpan(coRefPositions.breaks, forward,  coRefPositions.span, geneNames);
		int  span = geneNames.size();
		if(!TranscriptUtils.coronavirus && span==0 && q_value < ViralTranscriptAnalysisCmd2.fail_thresh1){
			return secondKey.toString();
		}
		String breakSt = coRefPositions.breaks.toString();
		//coRefPositions.breaks.adjustBreaks(annot);
		// need to group by start position if we annotating by break pos,e.g. so 5'mapping reads map together
		String secondKeySt = secondKey.toString();
		coRefPositions.breaks_hash.setSecondKey(secondKeySt);//, endPos);
		
		
		
		if(cluster_reads)  parent.all_clusters.matchCluster(coRefPositions, this.source_index, this.num_sources,  this.chrom_index, clusterID, strand); // this also clears current cluster
		else{
		clusterID[0] = chrom_+".NA";
		clusterID[1] = "NA";
		}
	
		String str = id+"\t"+clusterID[0]+"\t"+clusterID[1]+"\t"+source_index+"\t"+readLength+"\t"+start_read+"\t"+end_read+"\t"
		+type_nme+"\t"+chrom_+"\t"
		+startPos+"\t"+endPos+"\t"+(forward? "+":"-")+"\t"+coRefPositions.numIsoforms()+"\t"+(hasLeaderBreak ? 1:0)+"\t"
		+coRefPositions.getError(src_index)+"\t"+secondKeySt+"\t"+strand+"\t"+breakSt+"\t"+span_str+"\t"+geneNames.size()+"\t"+String.format("%5.3g", q_value).trim();
		parent.o.printRead(str);
		int num_exons =(int) Math.floor( (double)  coRefPositions.breaks.size()/2.0);

		boolean writeMSA = Outputs.doMSA!=null && Outputs.msa_sources !=null && includeInConsensus  && Outputs.msa_sources.containsKey(source_index);
		if(includeInConsensus && TranscriptUtils.coronavirus){
			int st1 = startPos; //position>0 ? position : startPos; // start after break
			inner: for(int i=annot.start.size()-1; i>=0; i--){
				if(st1 -10> annot.start.get(i)) break inner;
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
							parent.o.writeToCluster("annot_"+annot.genes.get(i),null, source_index, readSeq1, baseQ1,null, readSeq.getName(), strand);
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
			parent.o.writeToCluster(prefix+secondKeySt,"_"+clusterID[1]+"_", source_index, readSeq1, baseQ1, str, readSeq.getName(), strand);
		}
		return secondKeySt+" "+span_str;
	}
	
	

	private static String getString(Integer[] seq12) {
		return CigarHash2.getString(Arrays.asList(seq12));
	}
	

	public void addRefPositions(int position, boolean match) {
		//we add back in one here to convert it to a 1 based index
		coRefPositions.add(position+1, this.source_index, match);
	}


	static class CigarClusterRepo{
		final int num_sources;
		CigarClusterRepo(int num_sources){
			this.num_sources = num_sources;
		}
		Stack<CigarCluster> s = new Stack<CigarCluster>();
	
		public synchronized CigarCluster get(){
			if(s.size()==0){
				return new CigarCluster("-1",num_sources,'.');
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
	

	
	 

	
	
	


	/**
	 * Get the identity between a read sequence from a sam and a reference sequence
	 * 
	 * @param refSeq
	 * @param sam
	 * @return
	 */
	public  void identity1(Sequence refSeq, Sequence fivePrimeRefSeq, Sequence threePrimeRefSeq,   Sequence readSeq, SAMRecord sam, 
			int source_index, boolean cluster_reads,  String poolID, double qval) throws NumberFormatException{
		int readPos = 0;// start from 0
		int seqlen = refSeq.length();
		IdentityProfile1 profile = this;
		//Outputs output = profile.o;
		int refPos = sam.getAlignmentStart() - 1;// convert to 0-based index
		String id = sam.getReadName();
		String chrom = refSeq.getName();
		profile.updateSourceIndex(source_index);
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
				for (int i = 0; i < length && refPos + i < refSeq.length(); i++) {
				//	profile.baseDel[refPos + i] += 1;
					profile.addRefPositions(refPos + i, false);
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
				for (int i = 0; i < length && refPos + i < refSeq.length(); i++) {
					
					if (refSeq.getBase(refPos + i) == readSeq.getBase(readPos + i)){
				//		profile.match[refPos + i]++;
						profile.addRefPositions(refPos + i, true);
					}
					else{
					//	profile.mismatch[refPos + i]++;
						profile.addRefPositions(refPos + i, false);
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
				for (int i = 0; i < length && (refPos + i) < refSeq.length(); i++) {
					profile.addRefPositions(refPos+i, true);
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
				for (int i = 0; i < length && (refPos + i) < refSeq.length(); i++) {
					profile.addRefPositions(refPos+i, false);
				}
				//profile.addRefPositions(refPos);
			//	profile.mismatch[refPos] += length;
				break;
			default:
				throw new IllegalStateException("Case statement didn't deal with cigar op: " + e.getOperator());
			}// case
		} // for
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
			int polyAlen = TranscriptUtils.polyAlen(refSeq);
			if(TranscriptUtils.coronavirus && st_r > extra_threshold1 && attempt5rescue && st_r < 1000 && sam.getAlignmentStart()> 100 ){
				//good candidate for a missed 5' alignment
				try{
				 align_5prime = SWGAlignment.align(readSeq.subSequence(0, st_r), fivePrimeRefSeq);
				}catch(Exception exc){
					System.err.println("warning could not attemp 5'rescue "+ readSeq.subSequence(0, st_r));
				}
			}
			
			int seqlen1 = refSeq.length()-polyAlen;
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
				if(polyAlign.getIdentity() > 0.9 * polyAlign.getLength()  && polyAlign.getLength()>15){
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
			 
			String secondKey= profile.processRefPositions(sam, id, cluster_reads, 
					refSeq, source_index, readSeq,baseQ, st_r, end_r, strand, align_5prime, align_3prime,align_3primeRev, offset_3prime, polyAlen, poolID, qval);
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

				if(sam.getReadNegativeStrandFlag()) {
					leftseq =TranscriptUtils.revCompl(leftseq); 
					leftseq.setDesc(leftseq.getDesc()+" rev "+chrom);		
					baseQL = (new StringBuilder(baseQ)).reverse().toString();
				}
				else{
					leftseq.setDesc(leftseq.getDesc()+" fwd "+chrom);		
				}
			
				leftseq.setName(readSeq.getName());
				StringBuffer desc = new StringBuffer();
				desc.append("L "+secondKey+" "+st_r+" "+getString(seq11));
			
				if(reAlignExtra){
					Sequence refSeq1 = refSeq.length() < 30000  ? refSeq : 
							refSeq.subSequence(profile.startPos - 10000, profile.endPos+1000);
				align_5prime = SWGAlignment.align(leftseq, refSeq1);
				 TranscriptUtils.getStartEnd(align_5prime, seq1, seq2, 0, 0, sam.getReadNegativeStrandFlag());
				String secondKey1 =  profile.all_clusters.annot.nextUpstream(seq1[2], profile.chrom_index, forward)+";"+profile.all_clusters.annot.nextDownstream(seq1[3], profile.chrom_index,forward);
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
					String secondKey1 =  profile.all_clusters.annot.nextUpstream(seq1[2], profile.chrom_index,forward)+";"+profile.all_clusters.annot.nextDownstream(seq1[3], profile.chrom_index,forward);

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

	

	



	public void clear() {
		this.coRefPositions.clear();
		
	}







	

}