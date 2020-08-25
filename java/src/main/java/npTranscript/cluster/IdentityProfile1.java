package npTranscript.cluster;

import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.SortedSet;
import java.util.TreeSet;

import htsjdk.samtools.SAMRecord;
import japsa.bio.np.barcode.SWGAlignment;
import japsa.seq.Sequence;

/**
 * @author Lachlan Coin
 *
 */


public class IdentityProfile1 {
	

	final String[] type_nmes; 
	
	public static boolean annotByBreakPosition = true;
	
	public static int writeCoverageDepthThresh = 100;
	public static int[] writeIsoformDepthThresh = new int[] {10};
	public static int msaDepthThresh = 10;
	public static boolean includeStart = true;
	
	public IdentityProfile1(Sequence refSeq,
			Outputs o,
			String[] in_nmes,  int startThresh, int endThresh, boolean calcBreakpoints, Sequence chrom, int chrom_index) throws IOException {
	this.ref = chrom;
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
		all_clusters =new CigarClusters(refSeq,  num_sources);
	
		
		if(calcBreakpoints){
			System.err.println("calculating break point usage");
			bp = new BreakPoints(num_sources, seqlen);
		}else{
			bp  = null;
		}
		
	}
public void update(Annotation annot){
	this.all_clusters.update(annot);
}
	

	
	
	static double break_round = 10.0;
	static String NAstring = "NA";

	public static boolean subclusterBasedOnStEnd = false;
	public static int cleanUpMoreThan = 100000; //  distance between previous end and current start to remove cluster
	public static int clearUpThreshold = 100; //number of clusters before try to clear up 
	
	public int startPos, endPos;
	public int readSt, readEn; 
	
	SortedSet<String> geneNames = new TreeSet<String>();
	
	public String[] clusterID = new String[2];
	public String processRefPositions(SAMRecord sam, String id, boolean cluster_reads, Sequence refSeq, int src_index , Sequence readSeq, String baseQ, 
			int start_read, int end_read, char strand, SWGAlignment align5prime, SWGAlignment align3prime,
			SWGAlignment align3primeRev,
			int offset_3prime, int polyAlen, String pool
			) throws IOException, NumberFormatException{
		this.all_clusters.clearUpTo(sam.getAlignmentStart() -cleanUpMoreThan , o, refSeq, chrom_index, geneNames);
		startPos = sam.getAlignmentStart()+1; // transfer to one based
		endPos = sam.getAlignmentEnd()+1;
		readSt = start_read; readEn = end_read;
		boolean forward = !sam.getReadNegativeStrandFlag();
		//boolean hasSplice = false;
		int  readLength = readSeq.length();
		Annotation annot = this.all_clusters.annot;

		CigarHash2 breaks  = coRefPositions.breaks;
		int seqlen = refSeq.length();
	
		if(polyAlen>0){
			seqlen = seqlen-(polyAlen);
			annot.adjust3UTR(seqlen);
		}
		
		coRefPositions.end = endPos;
		coRefPositions.start = startPos;
		coRefPositions.forward = forward;
		breaks.add(coRefPositions.end);
		boolean includeInConsensus = true;
		if( align5prime!=null ){
			if(align5prime.getIdentity()>0.8 * Math.max(start_read,align5prime.getLength())){
				System.err.println("rescued 5' "+readSeq.getName());
			
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
			if(ident > 0.8 *Math.max(alignLen,diff)){// && read_st< 20){
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
			secondKey.append(annot.nextUpstream(startPos,chrom_index, forward)+";");
		}
		if(annotByBreakPosition){
			
			for(int i=1; i<breaks.size()-1; i+=2){
				int gap = breaks.get(i+1)-breaks.get(i);
				String upst = annot.nextUpstream(breaks.get(i), chrom_index,forward);
				secondKey.append(upst+",");
				secondKey.append(annot.nextDownstream(breaks.get(i+1), chrom_index,forward)+";");
				if(gap > TranscriptUtils.break_thresh){
					if(annot.isLeader(breaks.get(i))){
						if(bp!=null) this.bp.addBreakPoint(source_index, 0, breaks.get(i), breaks.get(i+1));
						hasLeaderBreak = true;
					}else{
						if(bp!=null)  this.bp.addBreakPoint(source_index, 1, breaks.get(i), breaks.get(i+1));

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
		
		String type_nme = annot.getTypeNme( startPos, endPos, forward); //coRefPositions.getTypeNme(seqlen);
		geneNames.clear();
		String span_str = annot.getSpan(coRefPositions.breaks, forward,  coRefPositions.span, geneNames);
		String breakSt = coRefPositions.breaks.toString();
		//coRefPositions.breaks.adjustBreaks(annot);
		// need to group by start position if we annotating by break pos,e.g. so 5'mapping reads map together
		String secondKeySt = secondKey.toString();
		coRefPositions.breaks_hash.setSecondKey(secondKeySt);
		
		
		
		if(cluster_reads)  this.all_clusters.matchCluster(coRefPositions, this.source_index, this.num_sources,  this.chrom_index, clusterID, strand); // this also clears current cluster
		else{
		clusterID[0] = chrom_+".NA";
		clusterID[1] = "NA";
		}
		int  span = geneNames.size();
		String str = id+"\t"+clusterID[0]+"\t"+clusterID[1]+"\t"+source_index+"\t"+readLength+"\t"+start_read+"\t"+end_read+"\t"
		+type_nme+"\t"+chrom_+"\t"
		+startPos+"\t"+endPos+"\t"+(forward? "+":"-")+"\t"+coRefPositions.numBreaks()+"\t"+(hasLeaderBreak ? 1:0)+"\t"
		+coRefPositions.getError(src_index)+"\t"+secondKeySt+"\t"+strand+"\t"+breakSt+"\t"+span_str+"\t"+geneNames.size();
		this.o.printRead(str);
		int num_exons =(int) Math.floor( (double)  coRefPositions.breaks.size()/2.0);

		this.o.printBed(coRefPositions.breaks, id,  strand, source_index, clusterID[0], clusterID[1], num_exons, span, span_str);
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
							this.o.writeToCluster("annot_"+annot.genes.get(i),null, source_index, readSeq1, baseQ1,null, readSeq.getName(), strand);
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
			this.o.writeToCluster(prefix+secondKeySt,"_"+clusterID[1]+"_", source_index, readSeq1, baseQ1, str, readSeq.getName(), strand);
		}
		return secondKeySt+" "+span_str;
	}
	
	

	

	

	public void addRefPositions(int position, boolean match) {
		//we add back in one here to convert it to a 1 based index
		coRefPositions.add(position+1, this.source_index, match);
	}

	
	final CigarCluster coRefPositions; // this is a temporary object used to store information for each read as it comes through 
	
	public CigarClusters all_clusters;
	final public BreakPoints bp;
	
	public void refresh(Sequence chr, int chrom_index){
		this.all_clusters.clear();
		this.ref = chr;
		this.chrom_ = chr.getName();
		this.chrom_index = chrom_index;
		if(bp!=null) bp.refresh(chr.length());
		this.genome =chr;
	}
	
	public int refBase, readBase;

	//private  PrintWriter readClusters;
	private Sequence genome;
	private final int num_sources;
	public int source_index=-1; //source index
	public void updateSourceIndex(int i) {
		this.source_index = i;
	}
	

	
	 

	public  Sequence ref;
	public String chrom_;
	public  int chrom_index;
	public final Outputs o;
	
	public void getConsensus() throws IOException {
			this.all_clusters.getConsensus( o, this.ref, this.chrom_index, this.geneNames);
			
	}	


	

	

	public void newRead(int source_index2) {
		this.updateSourceIndex(source_index2);
		this.coRefPositions.clear(source_index2);
	}







	public void printBreakPoints() throws IOException {
		if(bp!=null) this.bp.printBreakPoints(this.o, this.chrom_index);
		
	}

}