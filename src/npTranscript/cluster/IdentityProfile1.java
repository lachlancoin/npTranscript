package npTranscript.cluster;

import java.io.File;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.io.PrintWriter;
import java.util.Iterator;
import java.util.List;
import java.util.zip.GZIPOutputStream;

import org.apache.commons.math3.linear.OpenMapRealMatrix;
import org.apache.commons.math3.linear.SparseFieldMatrix;
import org.apache.commons.math3.linear.SparseRealMatrix;

import ch.systemsx.cisd.hdf5.HDF5Factory;
import ch.systemsx.cisd.hdf5.IHDF5SimpleWriter;
import japsa.seq.Sequence;
import japsa.seq.SequenceOutputStream;

/**
 * @author Lachlan Coin
 *
 */


public class IdentityProfile1 {
	final File outfile, outfile1, outfile2, outfile3, outfile4, outfile5, outfile6, outfile7, outfile8, outfile9;
	//OpenMapRealMatrix
	SparseFieldMatrix sm;
	
	
	
	
	public IdentityProfile1(Sequence refSeq, File resDir, List<Integer[]> positions, int num_sources, int genome_index,  boolean calculateCoExpression, double overlapThresh, int startThresh, int endThresh) throws IOException {
		//File readClusterFile = new File(outdir, "readclusters.txt.gz");
		if(positions.size()>0){
			int rowlen = refSeq.length();
			sm = new SparseFieldMatrix(new IField(), rowlen, rowlen);
		for(int i=0; i<positions.size(); i++){
			Integer[] pos = positions.get(i);
			sm.setEntry(pos[0], pos[1],  new IntegerField(pos[2]));
		}
		}
		
		//this.bin = (double) bin;
		this.num_sources = num_sources;
		this.coRefPositions = new CigarCluster(-1,num_sources);
//	 TranscriptUtils.round2 = 100.0/round;
		this.genome = refSeq;
		this.startThresh = startThresh; this.endThresh = endThresh;
		this.source_index = 0;		
		this.calculateCoExpression = calculateCoExpression;
		 outfile = new File(resDir,genome_index+ ".txt");
		 outfile1 = new File(resDir, genome_index+ "coref.txt");
		 outfile2 = new File(resDir, genome_index+"clusters.h5");
		 outfile3 = new File(resDir,genome_index+ "readToCluster.txt.gz");
		 outfile4 = new File(resDir,genome_index+ "exons.txt.gz");
		 outfile5 = new File(resDir,genome_index+ "clusters.fa.gz");
		 outfile6 = new File(resDir,genome_index+ "tree.txt.gz");
		 outfile7 = new File(resDir,genome_index+ "dist.txt.gz");
		 outfile8 = new File(resDir,genome_index+ "transcripts.txt.gz");
		 outfile9 = new File(resDir,genome_index+ "breakpoints.txt.gz");

		 readClusters = new PrintWriter(
					new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(outfile3))));
			//this.readClusters.println("readID,clusterID,index,source_index");//+clusterID+","+index+","+source_index);
		 
		readClusters.println("readID\tclusterId\tsource\tbreak_cluster\tlength\tstartPos\tendPos\tbreakStart\tbreakEnd\terrorRatio");
		
		/*
		readClipped = 0;
		numDel = 0;
		numIns = 0;
		match = new int[refSeq.length()];
		mismatch = new int[refSeq.length()];
		refClipped = new int[refSeq.length()];
		baseDel = new int[refSeq.length()];
		baseIns = new int[refSeq.length()];

		Arrays.fill(match, 0);
		Arrays.fill(mismatch, 0);
		Arrays.fill(baseDel, 0);
		Arrays.fill(baseIns, 0);
		Arrays.fill(refClipped, 0);
		 */
		refBase = 0;
		readBase = 0;
		// the number of bases from ref and read
		// following should be made more efficient
		int seqlen = refSeq.length();
		//Set<Integer> roundedPos = new HashSet<Integer>();
		//for (int i = 0; i < refSeq.length(); i++) {
		//	roundedPos.add(i);
		//}
	//	roundedPositions = roundedPos.toArray(new Integer[0]);
		//this.depth = new int[roundedPos.size()];
	
		all_clusters =new CigarClusters(overlapThresh);
		this.breakpoints = new SparseRealMatrix[this.num_sources];
		this.breakSt = new SparseVector[this.num_sources];
		this.breakEnd = new SparseVector[this.num_sources];
		for(int i=0; i<breakpoints.length; i++){
			this.breakSt[i] = new SparseVector();
			this.breakEnd[i] = new SparseVector();
			breakpoints[i] = new OpenMapRealMatrix(seqlen, seqlen);
		}
		
	}

	

	final int startThresh, endThresh;
	
	//final  double bin ;
	private final boolean calculateCoExpression;

//	final String[] nmes =  "5_3:5_no3:no5_3:no5_no3".split(":");

	//this after all postions in a read processed
	
	static double break_round = 10.0;
	public void processRefPositions(int startPos, int endPos, String id, boolean cluster_reads, int  readLength, int refLength, int src_index) {
		List<Integer> breaks  = coRefPositions.breaks;
		coRefPositions.end = endPos;
		coRefPositions.start = startPos;
		breaks.add(coRefPositions.end);
		int distToEnd = refLength - endPos;	
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
		if(maxg>100){
			prev_position = breaks.get(maxg_ind);
			position = breaks.get(maxg_ind+1);
			if(this.sm!=null){
				 coRefPositions.break_point_cluster = ((IntegerField)this.sm.getEntry(prev_position, position)).getVal();
				
			}
			coRefPositions.breakSt = prev_position;
			coRefPositions.breakEnd = position;
			this.breakpoints[this.source_index].addToEntry(prev_position, position, 1);
			this.breakSt[this.source_index].addToEntry(prev_position, 1);
			this.breakEnd[this.source_index].addToEntry(position,  1);
		}
		
		/*int index = 0;
		if (startPos < startThresh)
			index = distToEnd < endThresh ? 0 : 1;
		else
			index = distToEnd < endThresh ? 2 : 3;*/
	//	Iterator<Integer> it = coRefPositions.keys();
		coRefPositions.breaks.roundBreaks();
		
	//	coRefPositions.index = index;
		int clusterID;
		if(cluster_reads) clusterID = this.all_clusters.matchCluster(coRefPositions, this.source_index, this.num_sources); // this also clears current cluster
		else clusterID = -2;
		this.readClusters.println(id+"\t"+clusterID+"\t"+source_index+"\t"+coRefPositions.break_point_cluster+"\t"+readLength+"\t"
		+startPos+"\t"+endPos+"\t"+prev_position+"\t"+position+"\t"+coRefPositions.getError(src_index));
	}
	
	public void addRefPositions(int position, boolean match) {
		//we add back in one here to convert it to a 1 based index
		coRefPositions.add(position+1, this.source_index, match);
		
	}

	
	private final CigarCluster coRefPositions;
	private CigarClusters all_clusters;
	final public SparseRealMatrix[] breakpoints;
	final public SparseVector[] breakSt, breakEnd;
	//private Integer[] roundedPositions;// , corefSum;
	private int[] depth; // convenience matrix for depth
	//public int[] match, mismatch, refClipped, baseDel, baseIns;
	//public int numIns, numDel, readClipped,
	public int refBase, readBase;

	private final PrintWriter readClusters;
	private final Sequence genome;
	private final int num_sources;
	public int source_index; //source index
	public void updateSourceIndex(int i) {
		this.source_index = i;
	}
	
	/*public void print(PrintWriter pw, Sequence seq) {
		pw.println("pos,base,match,mismatch,refClipped,baseDel,baseIns");
		for (int i = 0; i < match.length; i++) {
			pw.println((i + 1) + "," + seq.charAt(i) + "," + match[i] + "," + mismatch[i] + "," + refClipped[i]
					+ "," + baseDel[i] + "," + baseIns[i]);
		}
		pw.flush();
	}*/

	/*public void printClusters(File outfile1) throws IOException {
		
			PrintWriter pw = new PrintWriter(
					new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(outfile1))));
			StringBuffer sb = new StringBuffer("pos,NA,"+getString("NA",num_sources,false));
			for (int i = 0; i < this.roundedPositions.length; i++) {
				sb.append(",");
				sb.append(roundedPositions[i] * round + 1);
			}
			pw.println(sb.toString());
			for (int i = 0; i < this.all_clusters.l.size(); i++) {
				this.all_clusters.l.get(i).print(pw);
			//	pw.println(this.all_clusters.l.get(i).summary(this.roundedPositions));
			}
			pw.close();
	}*/

	private void printBreakPoints(File outfile1)  throws IOException {
		// TODO Auto-generated method stub
		for(int i=0; i<this.breakpoints.length; i++){
		File outfile1_ = new File(outfile1.getParentFile(),
				outfile1.getName() + "." + i + ".gz");
		printMatrix(this.breakpoints[i],this.breakSt[i], this.breakEnd[i],  outfile1_);
		}
	}
	 void printMatrix(SparseRealMatrix cod, SparseVector breakSt2, SparseVector breakEnd2, File outfile1_) throws IOException{
		
		PrintWriter pw = new PrintWriter(
				new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(outfile1_))));
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
		pw.close();
	}

	
	
	public void getConsensus(Annotation annot) throws IOException {
			//int num_types = this.nmes.length;
			this.depth = new int[this.genome.length()];
			//for(int i=0; i<num_types; i++){
			
			PrintWriter exonsP =new PrintWriter( new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(outfile4))));
			
			SequenceOutputStream seqFasta =  new SequenceOutputStream(new GZIPOutputStream(new FileOutputStream(outfile5)));
			PrintWriter transcriptsP =  new PrintWriter( new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(outfile8))));
			IHDF5SimpleWriter clusterW = 
					 HDF5Factory.open(outfile2);
			//		new PrintWriter(
				//	new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(outfile2+"."+nmes[i]+".gz"))));
			//}
			this.all_clusters.getConsensus(annot, this.genome, exonsP, transcriptsP, seqFasta,clusterW, this.depth, this.num_sources);
		//	for(int i=0; i<num_types; i++){
			clusterW.close();
			exonsP.close();
			seqFasta.close();
			transcriptsP.close();
			
		
		
	}
/*	public void printTree() throws IOException{
		PrintWriter treeP =  new PrintWriter(
				new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(outfile6))));
		PrintWriter distP =  new PrintWriter(
				new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(outfile7))));
		System.err.println("calculating distance matrix..");
		DistanceMatrix dm = this.all_clusters.getDistanceMatrix( distP);
		System.err.println("..done");
		NeighborJoiningTree tree = new NeighborJoiningTree(dm);
		//treeP.print(tree.toString());
		NodeUtils.printNH(treeP, tree.getRoot(), true, false, 0, true);
		treeP.close();
		
		
		distP.close();
		
	}*/

	public void finalise()  throws IOException{
		this.readClusters.close();
		IdentityProfile1 pr1 = this;
		PrintWriter pw = new PrintWriter(new FileWriter(outfile));
	//	pr1.print(pw, genome);
		pw.close();
	//	pr1.printCoRef(outfile1);
		pr1.printBreakPoints(outfile9);
	//	pr1.printClusters(outfile2);
		System.out.println("========================= TOTAL ============================");
		
	}

	

	public void newRead(int source_index2) {
		this.coRefPositions.clear(source_index2);
	}

}