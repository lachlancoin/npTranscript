package npTranscript.cluster;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.zip.GZIPOutputStream;

import org.apache.commons.math3.linear.OpenMapRealMatrix;
import org.apache.commons.math3.linear.SparseRealMatrix;

import ch.systemsx.cisd.hdf5.HDF5Factory;
import ch.systemsx.cisd.hdf5.IHDF5SimpleWriter;
import japsa.seq.Sequence;
import japsa.seq.SequenceOutputStream;
import npTranscript.cluster.CigarCluster.Count;

/**
 * @author Lachlan Coin
 *
 */


public class IdentityProfile1 {
	
	public static class Outputs{
		public File transcripts_file;
		public File reads_file; 
		private final File outfile, outfile1, outfile2, outfile2_1, outfile4, outfile5, outfile6, outfile7, outfile10;
		//outfile9;
		SequenceOutputStream seqFasta =  null; // new SequenceOutputStream(new GZIPOutputStream(new FileOutputStream(outfile5)));
		//i
		
		 private PrintWriter transcriptsP,readClusters;
		final IHDF5SimpleWriter clusterW, altT;
		//PrintWriter readClusters;
		File resDir;
		
		String[] type_nmes;
		public void close() throws IOException{
			if(seqFasta!=null) seqFasta.close();
			
			transcriptsP.close();
			readClusters.close();
			clusterW.close();
			this.altT.close();
		}
		int genome_index=0;
	//	int currentIndex=0;
		
		/*public void updateChromIndex(int currentIndex2) throws IOException{
			this.genome_index = currentIndex2;
			this.newReadCluster(currentIndex2);
		}*/
		/*private void newReadCluster(int genome_index) throws IOException{
			 if(readClusters!=null) readClusters.close();
			
		}*/
		public Outputs(File resDir,  String[] type_nmes, boolean overwrite, int currentIndex) throws IOException{
			this.type_nmes = type_nmes;
			this.genome_index= currentIndex;
			 this.resDir = resDir;
			 int num_sources = type_nmes.length;
			 outfile = new File(resDir,genome_index+ ".txt");
			 outfile1 = new File(resDir, genome_index+ "coref.txt");
			 outfile2 = new File(resDir, genome_index+".clusters.h5");
			 outfile2_1 = new File(resDir, genome_index+"clusters1.h5");
			
			 outfile4 = new File(resDir,genome_index+ ".exons.txt.gz");
			 outfile5 = new File(resDir,genome_index+ ".clusters.fa.gz");
			 outfile6 = new File(resDir,genome_index+ ".tree.txt.gz");
			 outfile7 = new File(resDir,genome_index+ ".dist.txt.gz");
			 outfile10 = new File(resDir,genome_index+".isoforms.h5");
			
			 reads_file = new File(resDir,genome_index+ ".readToCluster.txt.gz");
			 readClusters = new PrintWriter(
					new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(reads_file))));
			 readClusters.println("readID\tclusterId\tsubID\tsource\tlength\ttype_nme\tchrom\tstartPos\tendPos\tbreakStart\tbreakEnd\terrorRatio\tupstream\tdownstream");
		
			 transcripts_file = new File(resDir,genome_index+ ".transcripts.txt.gz");
			//	newReadCluster(genome_index);

			 File[] f = new File[] {outfile1, outfile2,  outfile10, transcripts_file};

			 for(int i=0; i<f.length; i++){
				 if(overwrite && f[i].exists()){
					 System.err.println("deleteing "+f[i].getAbsolutePath());
					 f[i].delete();
				 }
				 else if(f[i].exists()){
					 throw new RuntimeException("will not overwrite file "+f[i].getAbsolutePath());
				 }
			 }
			 transcriptsP =  new PrintWriter( new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(transcripts_file))));
				String transcriptP_header = "ID\tchrom\tstart\tend\ttype_nme\tstartBreak\tendBreak\tisoforms\tleftGene\trightGene\ttotLen\tcountTotal\t"+TranscriptUtils.getString("count", num_sources,true)
				+"\t"+TranscriptUtils.getString("depth", num_sources, true)+"\t"+TranscriptUtils.getString("errors", num_sources, true)
				+"\t"+TranscriptUtils.getString("error_ratio", num_sources, true);
				StringBuffer nme_info = new StringBuffer();
				for(int i=0; i<type_nmes.length; i++) nme_info.append(type_nmes[i]+"\t");
				transcriptsP.println("#"+nme_info.toString());
				transcriptsP.println(transcriptP_header);
			
			
			
			List<String> str = new ArrayList<String>();
			str.add("subID"); //str.add("ẗype"); 
			for(int j=0; j<num_sources; j++){
				str.add("depth"+j);
			}
			
			altT = 
					 HDF5Factory.open(outfile10);
			altT.writeStringArray("header", str.toArray(new String[0]));
			str.set(0, "pos");
			for(int j=0; j<num_sources; j++){
				str.add("errors"+j);
			}
		
			clusterW = 
					 HDF5Factory.open(outfile2);
			clusterW.writeStringArray("header", str.toArray(new String[0]));
		
		//	clusterW.writeStringArray("header", str.toArray(new String[0]));
		}

		public PrintWriter getBreakPointPw( String chrom, int i)  throws IOException{
			
				File outfile1_ = new File(resDir,chrom+".breakpoints."+type_nmes[i]+"."+i +".txt.gz");
				PrintWriter pw = new PrintWriter(
					new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(outfile1_))));
			
			return pw;
		}
		
		public IHDF5SimpleWriter getH5Writer(){
			
			return clusterW;
		}

		
		public void printTranscript(String str){
			this.transcriptsP.println(str);
		}
		public void printRead(String string) {
			this.readClusters.println(string);
			
		}

		public void writeIntMatrix(String id, int[][] matr) {
			this.clusterW.writeIntMatrix(id, matr);
			
		}

		public void writeString(String id, Map<CigarHash2, Count> all_breaks, int num_sources) {
			int[][]str = new int[all_breaks.size()][];
			{
			Iterator<Entry<CigarHash2, Count>> it = all_breaks.entrySet().iterator();
			int maxl =-1;
			boolean all_equal  = true;
			int extra = num_sources+1;
			for(int i=0;  it.hasNext();i++){
				Entry<CigarHash2,Count> ch = it.next();
				CigarHash2 key = ch.getKey();
				int len = key.size();
				if(i==0) maxl = len;
				else if(len>maxl) {
					maxl = len;
					all_equal = false;
				}else if(len<maxl){
					len = maxl;
				}
				str[i] = new int[len+extra];
				str[i][0]  = ch.getValue().id();
				System.arraycopy(ch.getValue().count(), 0, str[i], 1, num_sources);
			//	str[i][1]  = ch.getValue().count();
				for(int j=0; j<key.size(); j++){
					str[i][j+extra] = key.get(j)*CigarHash2.round;
				}
				
			}
			if(!all_equal){
				for(int i=0; i<str.length; i++){
					if(str[i].length<maxl+extra){
						int[] str1 = new int[maxl+extra];
						System.arraycopy(str[i], 0, str1, 0, str[i].length);
						str[i] = str1;
					}
				}
			}
			}
			try{
//			this.altT.writeCompoundArray(id,str);
				this.altT.writeIntMatrix(id, str);
			}catch(Exception exc){
				exc.printStackTrace();
			}
			
		}

		
	}
	
	
	//OpenMapRealMatrix
	//SparseFieldMatrix sm;
	final String[] type_nmes; 
	
	public static boolean annotByBreakPosition = true;
	
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
	public void processRefPositions(int startPos, int endPos, String id, boolean cluster_reads, int  readLength, int refLength, int src_index , int seqlen) throws IOException, NumberFormatException{
		CigarHash2 breaks  = coRefPositions.breaks;
		Annotation annot = this.all_clusters.annot;
		coRefPositions.end = endPos;
		coRefPositions.start = startPos;
		breaks.addR(coRefPositions.end);
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
		this.o.printRead(id+"\t"+clusterID[0]+"\t"+clusterID[1]+"\t"+source_index+"\t"+readLength+"\t"
		+type_nme+"\t"+chrom+"\t"
		+startPos+"\t"+endPos+"\t"+prev_position+"\t"+position+"\t"+coRefPositions.getError(src_index)+"\t"+upstream+"\t"+downstream);
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