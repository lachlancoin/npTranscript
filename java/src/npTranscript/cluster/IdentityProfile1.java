package npTranscript.cluster;

import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import org.apache.commons.math3.linear.OpenMapRealMatrix;
import org.apache.commons.math3.linear.SparseRealMatrix;

import htsjdk.samtools.SAMRecord;
import japsa.bio.np.barcode.SWGAlignment;
import japsa.seq.Sequence;

/**
 * @author Lachlan Coin
 *
 */


public class IdentityProfile1 {
	
	//OpenMapRealMatrix
	//SparseFieldMatrix sm;
	final String[] type_nmes; 
	
	public static boolean annotByBreakPosition = true;
	
	public static int writeCoverageDepthThresh = 100;
	public static int writeIsoformDepthThresh = 10;
	public static int msaDepthThresh = 10;
	
	public IdentityProfile1(Sequence refSeq,
			Outputs o,
			String[] in_nmes,  int startThresh, int endThresh, Annotation annot, boolean calcBreakpoints, String chrom, int chrom_index) throws IOException {
	this.chrom = chrom;
	this.chrom_index = chrom_index;
	this.type_nmes = in_nmes;
		this.num_sources = in_nmes.length;
		this.coRefPositions = new CigarCluster("-1",num_sources);
		this.genome = refSeq;
		this.source_index = 0;		
		this.o  = o;
		refBase = 0;
		readBase = 0;
		int seqlen = refSeq.length();
		all_clusters =new CigarClusters(refSeq, annot, num_sources);
		this.breakpoints = new SparseRealMatrix[this.num_sources][2];
		this.breakSt = new SparseVector[this.num_sources][2];
		this.breakEnd = new SparseVector[this.num_sources][2];
		
		if(calcBreakpoints){
			System.err.println("calculating break point usage");
			for(int i=0; i<breakpoints.length; i++){
				this.breakSt[i] = new SparseVector[] {new SparseVector(),new SparseVector()};
				this.breakEnd[i] = new SparseVector[] {new SparseVector(),new SparseVector()};
				breakpoints[i] = new OpenMapRealMatrix[] {new OpenMapRealMatrix(seqlen, seqlen),new OpenMapRealMatrix(seqlen, seqlen)};
			
			}
		}
		
	}

	

	//final int startThresh, endThresh;
	
	//final  double bin ;
	//private final boolean calculateCoExpression;

//	final String[] nmes =  "5_3:5_no3:no5_3:no5_no3".split(":");

	//this after all postions in a read processed
	
	static double break_round = 10.0;
	static String NAstring = "NA";

	public static boolean subclusterBasedOnStEnd = false;
	
	public String[] clusterID = new String[2];
	public boolean processRefPositions(SAMRecord sam, String id, boolean cluster_reads, int  readLength, Sequence refSeq, int src_index , Sequence readSeq, String baseQ, 
			int start_read, int end_read, char strand, SWGAlignment align5prime
			) throws IOException, NumberFormatException{
		int startPos = sam.getAlignmentStart()+1; // transfer to one based
		int endPos = sam.getAlignmentEnd()+1;
		boolean hasSplice = false;
		Annotation annot = this.all_clusters.annot;

		CigarHash2 breaks  = coRefPositions.breaks;
		int seqlen = refSeq.length();
		char[] last10bp = refSeq.subSequence(seqlen-10, seqlen).charSequence();
		boolean polyA = true;
		for(int i=0; i<last10bp.length; i++){
			if(last10bp[i]!='A') polyA = false;
			
		}
		if(polyA){
			int i =1;
			while(refSeq.charAt(seqlen-i)=='A'){
				i++;
			}
			seqlen = seqlen-(i-1);
			annot.adjust3UTR(seqlen);
		}
		
		coRefPositions.end = endPos;
		coRefPositions.start = startPos;
		breaks.add(coRefPositions.end);
		
		if( align5prime!=null ){
			if(align5prime.getIdentity()>0.8 * Math.max(start_read,align5prime.getLength())){
				int newStartPos = align5prime.getStart2() + 1; // transfer to 1-based
				int newBreakPos = newStartPos + align5prime.getSequence2().length -  align5prime.getGaps2();
				if(newBreakPos < startPos){
					startPos = newStartPos;
					coRefPositions.start = newStartPos;
					coRefPositions.breaks.add(0, newBreakPos);
					coRefPositions.breaks.add(0, newStartPos);
				}
			}
		}
		//String roundStartP = annotByBreakPosition ? TranscriptUtils.round(startPos, CigarHash2.round)+"" : 	"";
		StringBuffer secondKey =new StringBuffer();
		//String upstream, upstream2, downstream, downstream2;
		int maxg = 0;
		int maxg_ind = -1;
		int maxg2=0;
		int maxg_ind2 =-1; 
		boolean hasLeaderBreak = breaks.size()>1 &&  annot.isLeader(breaks.get(1));
		secondKey.append(annot.nextUpstream(startPos,chrom_index)+";");
		if(annotByBreakPosition){
			
			for(int i=1; i<breaks.size()-1; i+=2){
				int gap = breaks.get(i+1)-breaks.get(i);
				String upst = annot.nextUpstream(breaks.get(i), chrom_index);
				secondKey.append(upst+",");
				secondKey.append(annot.nextDownstream(breaks.get(i+1), chrom_index)+";");
				if(gap > TranscriptUtils.break_thresh){
					if(annot.isLeader(breaks.get(i))){
						this.addBreakPoint(source_index, 0, breaks.get(i), breaks.get(i+1));
						hasLeaderBreak = true;
					}else{
						this.addBreakPoint(source_index, 1, breaks.get(i), breaks.get(i+1));

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
		secondKey.append(annot.nextDownstream(breaks.get(breaks.size()-1), chrom_index));
		
		
		
		String type_nme = coRefPositions.getTypeNme(seqlen);
		String breakSt = coRefPositions.breaks.toString();
		//coRefPositions.breaks.adjustBreaks(annot);
		// need to group by start position if we annotating by break pos,e.g. so 5'mapping reads map together
		String secondKeySt = secondKey.toString();
		coRefPositions.breaks_hash.setSecondKey(secondKeySt);
		
		if(cluster_reads)  this.all_clusters.matchCluster(coRefPositions, this.source_index, this.num_sources,  this.chrom_index, clusterID); // this also clears current cluster
		else{
		clusterID[0] = chrom+".NA";
		clusterID[1] = "NA";
		}
		
	//	System.err.println(id);
	//	String br_cluster_str = "";//sm==null ? "": coRefPositions.break_point_cluster+"\t";
		String str = id+"\t"+clusterID[0]+"\t"+clusterID[1]+"\t"+source_index+"\t"+readLength+"\t"+start_read+"\t"+end_read+"\t"
		+type_nme+"\t"+chrom+"\t"
		+startPos+"\t"+endPos+"\t"+coRefPositions.numBreaks()+"\t"+(hasLeaderBreak ? 1:0)+"\t"
		+coRefPositions.getError(src_index)+"\t"+secondKeySt+"\t"+strand+"\t"+breakSt;
		this.o.printRead(str);
		if(Outputs.doMSA!=null && (seqlen - endPos)<100){
			int st1 = startPos; //position>0 ? position : startPos; // start after break
			inner: for(int i=annot.start.size()-1; i>=0; i--){
			//	String gene = annot.genes.get(i);
				if(st1 -10> annot.start.get(i)) break inner;
				if(endPos > annot.end.get(i)){
					boolean include = false;
					inner1: for(int k=0; k<breaks.size(); k+=2){
						int st = breaks.get(k);
						int end = breaks.get(k+1);
						if(st-5 < annot.start.get(i) && end+5 > annot.end.get(i)){
							include=true;
							break inner1;
						}
					}
					if(include){
						int start_ref = annot.start.get(i)-1; // need to put back in zero coords
						int end_ref = annot.end.get(i)-1;// need to put back in zero coords
						int start_read1 =  sam.getReadPositionAtReferencePosition(start_ref, true);
						int end_read1 =  sam.getReadPositionAtReferencePosition(end_ref, true);
						double diff = end_ref - start_ref;
						double diff1 = end_read1 - start_read1;
						if(diff1>0 && diff1> 0.8 *diff){
							//int end_ref = sam.getReferencePositionAtReadPosition(end_read);
							//System.err.println(start_read1+","+end_read1+","+readSeq.length());
							//if(end_read1<start_read1){
							//	throw new RuntimeException("!!");
							//}
							Sequence readSeq1 = readSeq.subSequence(start_read1, end_read1);
							String baseQ1 = baseQ.length()<=1 ? baseQ : baseQ.substring(start_read1, end_read1);
							readSeq1.setDesc(chrom_index+";"+start_ref+","+end_ref+";"+start_read1+","+end_read1+";"+(end_read1-start_read1));
							this.o.writeToCluster("ORF_"+annot.genes.get(i),null, source_index, readSeq1, baseQ1,null, readSeq.getName(), strand);
						}
					}
				}
			}
		//	int end1 = endPos;
		}
		
		if(Outputs.doMSA!=null && Outputs.doMSA.contains(type_nme)) {
			Sequence readSeq1 = readSeq.subSequence(start_read, end_read);
			String baseQ1 = baseQ.length()<=1 ? baseQ :baseQ.substring(start_read, end_read);
		//	int end_ref = sam.getReferencePositionAtReadPosition(end_read);
		//	int start_ref = sam.getReferencePositionAtReadPosition(start_read);
			//readSeq1.setDesc("st="+start_read+";end="+end_read+";len="+(end_read-start_read));
			List<Integer>read_breaks = new ArrayList<Integer>();
			for(int i=0; i<breaks.size(); i++){
				read_breaks.add(sam.getReadPositionAtReferencePosition(breaks.get(i), true));
			}
			readSeq1.setDesc(chrom_index+";"+breaks.toString()+";"+CigarHash2.getString(read_breaks)+";"+(end_read-start_read));
			this.o.writeToCluster(secondKeySt,clusterID[1], source_index, readSeq1, baseQ1, str, readSeq.getName(), strand);
		}
		return hasSplice;
	}
	
	private void addBreakPoint(int source_index,int i, int prev_position, int position) {
		if(breakpoints[source_index][i]!=null){
			this.breakpoints[this.source_index][i].addToEntry(prev_position, position, 1);
			this.breakSt[this.source_index][i].addToEntry(prev_position, 1);
			this.breakEnd[this.source_index][i].addToEntry(position,  1);
		}
		
	}

	public void addRefPositions(int position, boolean match) {
		//we add back in one here to convert it to a 1 based index
		coRefPositions.add(position+1, this.source_index, match);
	}

	
	private final CigarCluster coRefPositions;
	private CigarClusters all_clusters;
	final public SparseRealMatrix[][] breakpoints;
	final public SparseVector[][] breakSt, breakEnd ;
	
	public int refBase, readBase;

	//private  PrintWriter readClusters;
	private final Sequence genome;
	private final int num_sources;
	public int source_index=-1; //source index
	public void updateSourceIndex(int i) {
		this.source_index = i;
	}
	

	public void printBreakPoints()  throws IOException {
		for(int i=0; i<this.breakpoints.length; i++){
			if(breakpoints[i]!=null){
				for(int j=0 ; j<breakpoints[i].length; j++)
				{
				PrintWriter pw = o.getBreakPointPw(chrom_index+"",i, j);
				printMatrix(this.breakpoints[i][j],this.breakSt[i][j], this.breakEnd[i][j],  pw);
				pw.close();
				}
			}
		}
	}
	 void printMatrix(SparseRealMatrix cod, SparseVector breakSt2, SparseVector breakEnd2, PrintWriter pw) throws IOException{
		
	
		StringBuffer secondLine = new StringBuffer();
		List<Integer> rows  = breakSt2.keys();
		List<Integer> cols =  breakEnd2.keys();
		for(Iterator<Integer> it = cols.iterator(); it.hasNext();){
			pw.print(",");
			Integer val = it.next();
			pw.print(val);
			secondLine.append(",");
			secondLine.append(breakEnd2.get(val));
		}
		pw.println();
		pw.println(secondLine.toString());
		//Set<Integer> cols = breakEnd2
		for (Iterator<Integer> it =rows.iterator(); it.hasNext();) {
			Integer row = it.next(); // nonZeroRows.get(i);
			pw.print(row);
			pw.print(",");
			pw.print(breakSt2.get(row));
			for(Iterator<Integer> it1 = cols.iterator(); it1.hasNext();){
				Integer col = it1.next();// nonZeroRows.get(j);
				int val =  (int) cod.getEntry(row, col);// : 0;
				pw.print(",");
				pw.print(val);
			}
			pw.println();
		}
		
	}

	public final String chrom;
	public final int chrom_index;
	public final Outputs o;
	
	public void getConsensus() throws IOException {
			this.all_clusters.getConsensus( o, this.chrom, this.chrom_index);
	}	


	

	

	public void newRead(int source_index2) {
		this.updateSourceIndex(source_index2);
		this.coRefPositions.clear(source_index2);
	}

}