package npTranscript.cluster;

import java.io.IOException;
import java.io.PrintWriter;
import java.util.Iterator;
import java.util.List;

import org.apache.commons.math3.linear.OpenMapRealMatrix;
import org.apache.commons.math3.linear.SparseRealMatrix;

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
		this.breakpoints = new SparseRealMatrix[this.num_sources];
		this.breakSt = new SparseVector[this.num_sources];
		this.breakEnd = new SparseVector[this.num_sources];
		if(calcBreakpoints){
			System.err.println("calculating break point usage");
			for(int i=0; i<breakpoints.length; i++){
				this.breakSt[i] = new SparseVector();
				this.breakEnd[i] = new SparseVector();
				breakpoints[i] = new OpenMapRealMatrix(seqlen, seqlen);
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
	
	String[] clusterID = new String[2];
	public void processRefPositions(int startPos, int endPos, String id, boolean cluster_reads, int  readLength, int refLength, int src_index , Sequence readSeq) throws IOException, NumberFormatException{
		CigarHash2 breaks  = coRefPositions.breaks;
		int seqlen = refLength;
		Annotation annot = this.all_clusters.annot;
		coRefPositions.end = endPos;
		coRefPositions.start = startPos;
		breaks.add(coRefPositions.end);
		//int distToEnd = refLength - endPos;	
		int maxg = 0;
		int maxg_ind = -1;
		for(int i=1; i<breaks.size()-1; i+=2){
			int gap = breaks.get(i+1)-breaks.get(i);
			if(gap>maxg){
				maxg = gap;
				maxg_ind = i;
			}
		}
		int prev_position = -1;
		int position = -1;
		String upstream = null;
		String downstream = null;
		if(maxg>100){
			prev_position = breaks.get(maxg_ind);
			position = breaks.get(maxg_ind+1);
			if(annotByBreakPosition ){
				if(annot!=null){
					downstream = annot.nextDownstream(position);
					upstream = annot.nextUpstream(prev_position);
				}
				if(upstream==null) upstream =  this.chrom_index+"."+TranscriptUtils.round(prev_position, CigarHash2.round);
				if(downstream==null) downstream = this.chrom_index+"."+ TranscriptUtils.round(position, CigarHash2.round);
			}
			/*if(this.sm!=null){
				 coRefPositions.break_point_cluster = ((IntegerField)this.sm.getEntry(prev_position, position)).getVal();
				
			}*/
			coRefPositions.breakSt = prev_position;
			coRefPositions.breakEnd = position;
			if(breakpoints[source_index]!=null){
				this.breakpoints[this.source_index].addToEntry(prev_position, position, 1);
				this.breakSt[this.source_index].addToEntry(prev_position, 1);
				this.breakEnd[this.source_index].addToEntry(position,  1);
			}
		}
		if(!annotByBreakPosition ){
			if(annot!=null){
				downstream = annot.nextDownstream(endPos);
				upstream = annot.nextUpstream(startPos);
			}
			if(upstream==null) upstream =  this.chrom_index+"."+TranscriptUtils.round(startPos, CigarHash2.round);
			if(downstream==null) downstream =  this.chrom_index+"."+TranscriptUtils.round(endPos, CigarHash2.round);

		}
		
		String type_nme = coRefPositions.getTypeNme(seqlen);
		//String breakSt1 = coRefPositions.breaks.toString();
		//coRefPositions.breaks.adjustBreaks(annot);
		coRefPositions.breaks_hash.setSecondKey(type_nme+"."+upstream+"."+downstream);
		
		if(cluster_reads)  this.all_clusters.matchCluster(coRefPositions, this.source_index, this.num_sources,  this.chrom_index, clusterID); // this also clears current cluster
		else{
		clusterID[0] = chrom+".NA";
		clusterID[1] = "NA";
		}
	//	System.err.println(id);
	//	String br_cluster_str = "";//sm==null ? "": coRefPositions.break_point_cluster+"\t";
		String str = id+"\t"+clusterID[0]+"\t"+clusterID[1]+"\t"+source_index+"\t"+readLength+"\t"
		+type_nme+"\t"+chrom+"\t"
		+startPos+"\t"+endPos+"\t"+prev_position+"\t"+position+"\t"+coRefPositions.getError(src_index)+"\t"+upstream+"\t"+downstream;
		this.o.printRead(str);
		this.o.writeToCluster(clusterID[0]+","+chrom+","+upstream+","+downstream+","+type_nme, source_index, readSeq, str);
	}
	
	public void addRefPositions(int position, boolean match) {
		//we add back in one here to convert it to a 1 based index
		coRefPositions.add(position+1, this.source_index, match);
	}

	
	private final CigarCluster coRefPositions;
	private CigarClusters all_clusters;
	final public SparseRealMatrix[] breakpoints;;
	final public SparseVector[] breakSt, breakEnd;
	
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
				PrintWriter pw = o.getBreakPointPw(chrom_index+"",i);
				printMatrix(this.breakpoints[i],this.breakSt[i], this.breakEnd[i],  pw);
				pw.close();
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