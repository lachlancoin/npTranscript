package npTranscript.cluster;

import java.awt.Color;
import java.io.BufferedReader;
import java.io.EOFException;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.OutputStream;
import java.io.OutputStreamWriter;
import java.io.PrintStream;
import java.io.PrintWriter;
import java.net.URI;
import java.net.http.HttpClient;
import java.net.http.HttpClient.Version;
import java.net.http.HttpRequest;
import java.net.http.HttpRequest.BodyPublishers;
import java.net.http.HttpResponse;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.SortedSet;
import java.util.TreeSet;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.ThreadPoolExecutor;
import java.util.stream.Collectors;
import java.util.stream.Stream;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;

import org.apache.commons.math3.linear.SparseRealMatrix;

import com.google.gson.Gson;

import ch.systemsx.cisd.hdf5.IHDF5SimpleWriter;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.fastq.FastqRecord;
import htsjdk.samtools.fastq.FastqWriter;
import htsjdk.samtools.fastq.FastqWriterFactory;
import japsa.seq.Sequence;
import npTranscript.NW.PolyAT;
import npTranscript.cluster.CigarCluster.Count;
import npTranscript.cluster.MultiSAMRecord.StringBuffers;
import npTranscript.run.Barcodes;
import npTranscript.run.CompressDir;
import npTranscript.run.SequenceOutputStream1;
import npTranscript.run.ViralTranscriptAnalysisCmd2;

public class Outputs{
	public static String readsOutputFile=null;
	public static boolean gzip =false;
	public static String url="http://0.0.0.0:81";
	public static PrintStream outputstream;  // this is the main output stream
	public static PrintStream joinOut,overlapOut,noGap,allOut, spliceOut, fiveOut, threeOut, remOut;
	public static String format;
	public static String sampleName;
	public static List<String> annotation_mode;
	
	//public static  ExecutorService executor ;
	
	public static final FastqWriterFactory factory = new FastqWriterFactory();
	{
		factory.setUseAsyncIo(true);
	}
	
	public static boolean overwrite = true;

	public static final SAMFileWriterFactory factB = new SAMFileWriterFactory();
	//public static int gffThreshGene = 10;
	
	public static boolean writeIsoforms = false;

	
	public static File library = new File("./");
	public static final ExecutorService writeCompressDirsExecutor  = Executors.newSingleThreadExecutor();
	//public static final ExecutorService fastQwriter = Executors.newSingleThreadExecutor();
	public static final ExecutorService h5writer = Executors.newSingleThreadExecutor();
	
	public static void shutDownExecutor() {
		System.err.println("shutting down output executors");
		if(writeCompressDirsExecutor!=null){
			waitOnThreads(writeCompressDirsExecutor,100);
			writeCompressDirsExecutor.shutdown();
		}
	/*	if(fastQwriter!=null){
			waitOnThreads(fastQwriter,100);

			Outputs.fastQwriter.shutdown();
		}*/
		if(h5writer!=null){
			waitOnThreads(h5writer,100);

			Outputs.h5writer.shutdown();
		}
		
	}
	
	
	/* class FOutp{
		OutputStreamWriter os;
		FastqWriter fastq ;
		File f; 
		FOutp(String nme	) throws IOException{
			this(nme, false);
		}
		FOutp(String nme, boolean fq) throws IOException{
			boolean gz = gzipFasta;
		//f = 
			if(fq) {
				f = new File(resDir,genome_index+nme+".fastq");
				fastq = factory.newWriter(f);
			}
			else{
				f = new File(resDir,nme+".fasta"+(gz ? ".gz": ""));
				OutputStream os1 = new FileOutputStream(f);
				if(gz) os1 = new GZIPOutputStream(os1);
				os =new OutputStreamWriter(os1);
			}
		}
		public void close() throws IOException {
			if(os!=null) this.os.close();
			else fastq.close();
			if(f.length()==0) f.deleteOnExit();
		}
	}*/
	public static List<String> doMSA = null;//"5_3";
	public static Map<Integer,Integer>msa_sources = new HashMap<Integer, Integer>();
//	public static boolean mergeSourceClusters = true;
	public static boolean gzipFasta = false;
	public static boolean keepAlignment = true;
	public static boolean keepinputFasta = true;
	//public static boolean writeBed=false;
	
	public static int minClusterEntries = 5;
	public static Collection numExonsMSA = Arrays.asList(new Integer[0]); // numBreaks for MSA 
	public static boolean calcBreaks;
	public static boolean writeH5=false;
	
		
		//public File reads_file; 
		
	//	private final File  outfile2,  outfile10,  outfile12;
		//private final File[] gff_output;
		//outfile9;
		//private final FOutp[] leftover_l, polyA;//, leftover_r, fusion_l, fusion_r;
	
	//	final int seqlen;
		//public  PrintWriter readClusters;//, plusMinus;
		//public  PrintWriter unmapped;//, plusMinus;

		// IHDF5SimpleWriter clusterW = null;
		 //IHDF5Writer altT = null;
		 //IHDF5Writer breakPW = null;
		static File resDir;
		CompressDir[] clusters;
		
		
		
		String type_nmes;
		public void close() throws IOException{
			System.err.println();
			this.all_res.values().stream().forEach(t ->t.values().forEach(t1 -> t1.post()) );
			Outputs.outputstream.close();
			Outputs.ps.stream().forEach(t->t.close());
			CompressDir cd = new CompressDir(resDir, false, gzip);
			cd.writeAll();
			cd.close();
			if(writeCompressDirsExecutor!=null){
				Outputs.waitOnThreads(writeCompressDirsExecutor,100);
				
			}
			/*if(fastQwriter!=null){
				Outputs.waitOnThreads(fastQwriter, 100);
			}*/
			if(h5writer!=null){
				Outputs.waitOnThreads(h5writer, 100);
			}
			//readClusters.close();
			//unmapped.close();
			/*if(clusterW!=null){
				writeAllDepths();
				clusterW.close();
			}*/
			//if(altT!=null) this.altT.close();
			//this.clusters.close();
			for(int i=0; i<clusters.length; i++){
				if(clusters[i]!=null) this.clusters[i].run(Outputs.minClusterEntries *2, writeCompressDirsExecutor);
			}
			
			
		

			
		}
		
		
			
		
		boolean writeDirectToZip = false;
		
		String genome_index;
		 int col_inds, col_inds_depth;  // this is the col_inds for writing count information in h5 library
		 int new_max_cols, new_max_cols_depth;
		private String last_id;
		public Outputs( String type_nmes, boolean isoforms, boolean cluster_depth, Map<String, Object> expt) throws IOException{
			this.type_nmes = type_nmes;
			this.expt = expt;
			this.register();
			genome_index = "0.";
		//	 this.resDir = resDir;
			 int num_sources = 1;//type_nmes.length;
			// outfile2 = new File(resDir, genome_index+"clusters.h5");
			 //outfile10 = new File(resDir,genome_index+"isoforms.h5");
			//outfile12 = new File(library,genome_index+"breakpoints.h5");
			List<Integer> vals = new ArrayList<Integer>(new HashSet<Integer> (Outputs.msa_sources.values()));
			Collections.sort(vals);
			 if(doMSA!=null && ( Outputs.msa_sources.size()==0)){
				 clusters = new CompressDir[] {new CompressDir(new File(resDir,  genome_index+"clusters"), true, gzip)};
			 }else if(doMSA!=null &&  ( Outputs.msa_sources.size()>0)){
				 clusters =  new CompressDir[vals.size()];
				 for(int i=0; i<clusters.length; i++){
					 String nmei =  genome_index+vals.get(i)+".";
					 clusters[i] = new CompressDir(new File(resDir, nmei+"clusters"), true, gzip);
				 }
			 }else{
				clusters = new CompressDir[1];
			 }
		
		//	this.right=  new SequenceOutputStream((new FileOutputStream(new File(resDir,genome_index+".right" ))));
		//	 reads_file = readsOutputFile==null ? null : new File(readsOutputFile);
			// if(!overwrite && reads_file.exists()){
			//	last_id = transferToNewFile(reads_file);//this makes sure there are no partial lines at the end 
			 //}else{
				 last_id = null;
					// TODO Auto-generated method stub
			// }
			// OutputStream os = new FileOutputStream(reads_file, !overwrite);
			// boolean gzip = true;
			 //if(reads_file.getName().endsWith(".gz")) os = new GZIPOutputStream(os);
			// readClusters = new PrintWriter(	new OutputStreamWriter(os));
				
			 
			// OutputStream os1 = null;//new FileOutputStream(reads_file+".unmapped.gz", !overwrite);
				// boolean gzip = true;
//				os1 = new GZIPOutputStream(os1);
	//			 unmapped = new PrintWriter(
	//					new OutputStreamWriter(os1));
			 
			
//			 readID  clusterId       subID   source  length  start_read      end_read   
			 //type_nme        chrom   startPos        endPos  breakStart      breakEnd        errorRatio
			 //upstream        downstream      strand  breaks

/*if(IdentityProfile1.trainStrand){
	header = header+ "\tleft10bp\tright10bp\tqual_left10bp\tqual_right10bp\tcorrect_strand";
}*/
			 String str1 = 	 ViralTranscriptAnalysisCmd2.barcodes == null ? "": Barcodes.getHeader();
;
//if(last_id==null){
if(!Outputs.format.equals("json")) {
			 Outputs.outputstream.println(IdentityProfile1.header+"\t"+PolyAT.getHeader()+"\t"+str1); //\tbreakStart\tbreakEnd\tbreakStart2\tbreakEnd2\tstrand\tbreaks");
}
			
				
			List<String> str = new ArrayList<String>();
			//str.add("subID"); //str.add("ẗype"); 
			//for(int j=0; j<num_sources; j++){
				str.add(this.type_nmes);
			//}
		//	this.col_inds = new int[this.type_nmes.length];
		/*	if(Outputs.writeH5){
				if(overwrite) outfile10.delete(); // remove existing h5 (temporary)
				altT = HDF5Factory.open(outfile10);
			}*/
/*			if(Outputs.calcBreaks && Outputs.writeH5){
				this.breakPW =  HDF5Factory.open(outfile12);

			}*/
			/*
			if(altT!=null) {
				writeStringArray(altT, "header", str.toArray(new String[0]),this.col_inds,0);
			}*/
			this.new_max_cols = col_inds +1;//Arrays.stream(col_inds) .max() .getAsInt()+1;
				     
				     
					
			
			str .clear();
			str.add("pos"); //str.add("ẗype"); 
		//	str.add("base");
//			for(int j=0; j<num_sources; j++){
				str.add(type_nmes);
	//		}
		/*	if(cluster_depth && writeH5){
				if(overwrite) outfile2.delete();
			clusterW = 	 HDF5Factory.open(outfile2);
			col_inds_depth =new int[num_sources];
			writeStringArray(clusterW, "header", str.toArray(new String[0]),col_inds_depth,1);
			this.new_max_cols_depth = Arrays.stream(col_inds_depth) .max() .getAsInt()+1;
			}*/
			
		//	clusterW.writeStringArray("header", str.toArray(new String[0]));
		}

		
		
		



		private String transferToNewFile(File reads_file) throws IOException {
			File newF = new File(reads_file.getParentFile(),"tmp."+System.currentTimeMillis()+".gz");
			OutputStream os = new GZIPOutputStream(new FileOutputStream(newF));
			PrintWriter pw = new PrintWriter(new OutputStreamWriter(os));
			InputStream is = new FileInputStream(reads_file);
			if(reads_file.getName().endsWith(".gz")) is = new GZIPInputStream(is);
			BufferedReader br = new BufferedReader(new InputStreamReader(is));
			String st = br.readLine();
			pw.println(st);
			String[] head = st.split("\t");
			int len = head.length;
			String last=null;
			while(st!=null){
				try{
					st = br.readLine();
					String[] str = st.split("\t");
					if(str.length!=head.length){
						break;
					}
					pw.println(st);
					last=str[0];
				}catch(EOFException exc){
					st = null;
					break;
				}
			
			}
			pw.close();
			br.close();
			String path = reads_file.getAbsolutePath();
			reads_file.delete();
			Files.move(Paths.get(newF.getAbsolutePath()), Paths.get(path));
			System.out.println("last usable is  "+last);
			return last;
			
		}







		/*public IHDF5SimpleWriter getH5Writer(){
			
			return clusterW;
		}*/

		/*
		private static void writeStringArray(IHDF5SimpleWriter altT2, String string, String array, int inds,  int start) {
			boolean exists = altT2.exists("header");
			List<String> newh= new ArrayList<String>();
			if(exists){
				String[] existing =altT2.readStringArray("header");
				newh.addAll(Arrays.asList(existing));
			}
			else{
				newh.addAll(Arrays.asList(array));
			}
			for(int i=start; i<array.length; i++){
			if(newh.contains(array[i])){
				if(inds!=null) inds[i-start] = newh.indexOf(array[i]);
			}else {
						
					if(inds!=null)		inds[i-start] = newh.size();
							newh.add(array[i]);
						}
				}
			altT2.writeStringArray(string, newh.toArray(new String[0]));
		}*/
		/*public synchronized void printTranscriptAlt(CigarCluster cc){
			if(altT!=null){
			String id2 = "trans/"+cc.breaks_hash.secondKey;
			int[]  obj;
			if(altT.exists(id2)){
				 obj = expand(altT.readIntArray(id2),this.new_max_cols);
			}else{
				obj= new int[this.new_max_cols];
				//obj = new HDFObjT(cc.breaks_hash.secondKey,cnts);
				//cigarClusters.infoString(cc, obj, geneNames, chrom);
			}
			for(int i=0; i<this.col_inds.length; i++){
				obj[col_inds[i]] = cc.readCount[i];
			}
			altT.writeIntArray(id2, obj);
			}
		}*/
		
		public synchronized void printRead(String string) {
			//System.err.println(string);
			if(!Outputs.format.equals("json")) {
					Outputs.outputstream.println(string);
			}
//			this.readClusters.println(string); this.readClusters.flush();
			
		}
		
		private PrintWriter getBreakPointPw( String chrom,  int j)  throws IOException{
			System.err.println("writing breakpoint files");
				File outfile1_ =  new File(resDir,chrom+".breakpoints."+type_nmes+"."+j +".txt.gz") ;
//						new File(resDir,chrom+".breakpoints."+type_nmes[i]+"."+i +".txt.gz");
				PrintWriter pw = new PrintWriter(
					new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(outfile1_, !overwrite))));
			
			return pw;
		}
		
		void printMatrix(SparseRealMatrix cod, SparseVector breakSt2, SparseVector breakEnd2,   int chrom_index, int j, String prefix) throws IOException{
			String id = prefix+"/"+this.type_nmes+"/"+j;
			//String id_scores = scores+"/"+this.chrom+"/"+this.type_nmes[i]+"/"+j;
//			PrintWriter pw = getBreakPointPw(chrom_index+"",i, j);
			StringBuffer secondLine = new StringBuffer();
			List<Integer> rows  = breakSt2.keys();
			List<Integer> cols =  breakEnd2.keys();
			int[][] mat = new int[rows.size()+2][cols.size()+2];
			for(int k=0; k<cols.size();k++){
				Integer val =  cols.get(k);
				mat[0][k+2] =val;
				mat[1][k+2] = breakEnd2.get(val);
			}
			//Set<Integer> cols = breakEnd2
			for (int m=0; m<rows.size(); m++) {
				Integer row = rows.get(m); // nonZeroRows.get(i);
				mat[m+2][0] = row;
				mat[m+2][1] = breakSt2.get(row);
				for(int k=0; k<cols.size();k++){
					Integer col = cols.get(k);
					int val =  (int) cod.getEntry(row, col);// : 0;
					mat[m+2][k+2] = val;
				}
				//pw.println();
			}
		//	this.breakPW.writeIntMatrix(id, mat);
			//pw.close();
			
		}
		
		static Color[] cols = new Color[] {Color.BLACK, Color.blue, Color.GREEN, Color.CYAN, Color.YELLOW, Color.red, Color.ORANGE, Color.MAGENTA, Color.pink};
		static int col_len = cols.length;
		static String[] col_str = new String[cols.length];

		//public static boolean all_isoforms = false;
		public static double isoThresh = 0.9;
		static{
			float[] f = new float[3];
			for(int i=0; i<col_str.length; i++){
				cols[i].getRGBColorComponents(f);
				
				col_str[i] = (int) Math.floor(255f*f[0])+","+(int) Math.floor(255f*f[1])+","+(int) Math.floor(255f*(f[2]));
			}
			//System.err.println(Arrays.asList(col_str));
		}
		
		
		
		public synchronized void writeToCluster(String ID, String subID,  int i, Sequence seq, String baseQ,  String name, char strand) throws IOException{
			CompressDir cluster = this.getCluster(i);
			//if(strand=='-'){
				//seq = TranscriptUtils.revCompl(seq);
				//baseQ = new StringBuilder(baseQ).reverse().toString();
			//}
		 	seq.setName(name);
			String entryname =  ID;
			if(subID!=null) entryname = entryname+"."+subID;
			if(cluster!=null){
				SequenceOutputStream1 out   = cluster.getSeqStream(entryname+".fa");
				//Outputs.writeFasta(out, seq);
				out.write(seq);
				//cluster.closeWriter(out);
			}
		}
		
		
		
		
			public synchronized void  writeFastq(FastqWriter writer, Sequence subseq,String baseQ,  boolean negStrand, int source_index)  throws IOException{
				if(writer==null) return;
						 writer.write(new FastqRecord(subseq.getName()+ " "+subseq.getDesc(), 
								 
								 new String( subseq.charSequence()), "", 
								 negStrand ? new StringBuilder(baseQ).reverse().toString() : baseQ));
			}

			
		Map<String, Maps> all_maps = new HashMap<String, Maps>();
		
		public synchronized void addDepthMap(String key, Maps map){
			if(map!=null){
			Maps m = all_maps.get(key);
			if(m==null) all_maps.put(key, new Maps(map));
			else {
				m.merge(map);
			//	all_maps.put(key, m);
			}
			}
		}
		
		/*public void writeAllDepths(){
			for(Iterator<String> keys = all_maps.keySet().iterator(); keys.hasNext();){
				String key = keys.next(); 
				Maps map = all_maps.get(key);
				writeDepthH5(key, map);
				keys.remove();	
			}
		}*/
		
		//	cc.breaks_hash.secondKey;
		/*
		public void writeDepthH5(String secondKey,Maps map) {
			int offset=1;
			int totalDepth = map.totalDepth();
			//Maps map = cc.map;
			if(clusterW!=null && totalDepth>=IdentityProfile1.writeCoverageDepthThresh){
				String  key = "depth/"+secondKey;
				
				

				//
				//cc.addZeros(cigarClusters.seqlen); 
				{
					List<Integer> keys = map.map.keys();

					int[][]matr = getMatr(key, offset, keys); //int[][] matr_err;
					map.getClusterDepth(matr, keys,  this.col_inds_depth, offset);
					clusterW.writeIntMatrix(key, matr);
				}
				if(CigarCluster.recordStartEnd){
					{
						String  keyStart = "depthStart/"+secondKey; // this for start end
						List<Integer> keys = map.mapStart.keys();
						int[][]matr = getMatr(keyStart, offset, keys); //int[][] matr_err;
						Maps.getClusterDepthStartEnd(matr, keys, offset,map.mapStart);
						clusterW.writeIntMatrix(keyStart, matr);
					}
					{
						String  keyEnd= "depthEnd/"+secondKey;
						List<Integer> keys = map.mapEnd.keys();
						int[][] matr = getMatr(keyEnd, offset, keys); //int[][] matr_err;
						Maps.getClusterDepthStartEnd(matr, keys, offset,map.mapEnd);
						clusterW.writeIntMatrix(keyEnd, matr);
					}
				}
				
			}
		}*/
		
	
		/*
		private int[][] getMatr( String key, int offset, List<Integer>keys) {
			int[][]matr;
			if(clusterW.exists(key)){
				int[][]  matr1 = clusterW.readIntMatrix(key);
				for(int i=0; i<matr1.length; i++){
					if(!keys.contains(matr1[i][0])){
						keys.add(matr1[i][0]);
					}
				}
				Collections.sort(keys);
				matr = new int[keys.size()][2*(new_max_cols_depth-offset)+offset];
				for(int i=0; i<matr1.length; i++){
					int ind = keys.indexOf(matr1[i][0]);
					System.arraycopy(matr1[i], 0, matr[ind], 0, matr1[i].length);
				}
				
			}else{
				matr = new int[keys.size()][2*(new_max_cols_depth-offset)+offset];
			}
			return matr;
		}

*/





		private synchronized CompressDir getCluster(int source){
			if(clusters==null) return null;
			return this.clusters[this.msa_sources.get(source)];
			//return mergeSourceClusters? this.clusters[0] : this.clusters[source];
		}
		private int getCount(Count val, int source){
			return val.count()[this.msa_sources.get(source)];
		//	return mergeSourceClusters ? Arrays.stream(val.count()).sum()  : val.count()[source];
		}

		public static String getString(int[] l,char sep){
			StringBuffer sb = new StringBuffer();
			for(int i=0; i<l.length; i++){
				if(i>0) sb.append(sep);
				sb.append(l[i]);
			}
			return sb.toString();
		}
		
		public static String getString(float[] l,char sep){
			StringBuffer sb = new StringBuffer();
			for(int i=0; i<l.length; i++){
				if(i>0) sb.append(sep);
				sb.append(l[i]);
			}
			return sb.toString();
		}
		public static class HDFObj  {
			public  String toString(){
				return getString(br,',')+"\t"+getString(cnts,',');//+"\t"+this.break_sum;
			}
			public HDFObj(String str) {
				String[] str_ = str.split("\t");
				//this.nme = str_[0];
				String[] str1 = str_[0].split(",");
				String[] str2 = str_[1].split(",");
				//this.break_sum=  Integer.parseInt(str_[2]);
				br = new float[str1.length];
				cnts = new int[str2.length];
				for(int i=0; i<str1.length; i++)br[i] = Float.parseFloat(str1[i]);
				for(int i=0; i<str2.length; i++)cnts[i] = Integer.parseInt(str2[i]);
			}
			public HDFObj(List<Integer> l1, int source_index, int new_max_cols, double divisor) {
				 cnts = new int[new_max_cols];
				cnts[source_index]=1;
					br = new float[l1.size()];
					for(int j=0; j<br.length; j++){
						br[j] = (float) (l1.get(j).doubleValue()/divisor);
					}
				//	this.break_sum = 1;
			}
			
			 float[]  br; int[] cnts;// int break_sum;
			/*@Override
			public int compareTo(Object o) {
				return nme.compareTo(((HDFObj)o).nme);
			}*/
			public void expand(int new_num_cols) {
				if(cnts.length < new_num_cols)
					cnts = Outputs.expand(cnts, new_num_cols);
			}
		}
		private static int[] expand(int[] cnts_old, int new_num_cols){
			int[] cnts = new int[new_num_cols];
			System.arraycopy(cnts_old, 0, cnts, 0, cnts_old.length);
			return cnts;
		}
		/*static class HDFObjT  {
			String key; int[] cnts;
			String chrom; int start; int end; String type_nme; String exon_count; int span; String genes;
		
			public HDFObjT(){
				
			}
			public HDFObjT(String key, int[] count) {
				this.key = key;
				cnts = count;
			}
			
			public void expand(int new_num_cols) {
				cnts = Outputs.expand(cnts, new_num_cols);
			}
		}*/
		
	//	String[] transcript_info  = new String[12];
		SortedSet<String> geneNames = new TreeSet<String>();
		
/*		public synchronized void writeTotalString(CigarClusters cigarClusters){
			if(altT==null) return;
			 String id2 = "counts/";//+cigarClusters.chr;
			int[]   obj;
			if(altT.exists(id2)){
				obj = expand( altT.readIntArray(id2),this.new_max_cols);
				
			}else{
				obj = new int[this.new_max_cols];
			
			}
			for(int i=0; i<this.col_inds.length; i++){
				obj[col_inds[i]] = cigarClusters.totalCounts[i];
			}
			altT.writeIntArray(id2, obj);
			
		}*/
		
		public static void main(String[] args) {
			/*
			
			int[] ints = {1, 2, 3, 4, 5};
			String[] strings = {"abc", "def", "ghi"};
Map<String,Map<String, int[]>> a = new HashMap<String, Map<String,int[]>>();
Map<String, int[]> ab=new HashMap<String, int[]>();
ab.put("aaab", new int[] {2,3});
a.put("aa", ab);
String str2 = gson.toJson(a);
			
			// Serialization
			String str1 = gson.toJson(ints);     // ==> [1,2,3,4,5]
			gson.toJson(strings);  // ==> ["abc", "def", "ghi"]
System.err.println(str2);
gson.fromJson(str1,  int[].class);
			// Deserialization
			int[] ints2 = gson.fromJson("[1,2,3,4,5]", int[].class);*/
		}

		
		
	//	Map<String,Map<String,Map<String,Map<String,Map<String,Object>>>>> all_res = 
	//			new HashMap<String,Map<String,Map<String,Map<String,Map<String,Object>>>>>();
		
		Map<String,Map<String,JSONOut>> all_res = new HashMap<String,Map<String, JSONOut>>(); //chrom strand
		
		class JSONOut{
			String chr;String strand;
			String first_read;
			//boolean transplice;
			boolean annotate;
			Set<String> reads = new TreeSet<String>();
			Map<String,Map<String,Map<String,Object>>> all_res =  //round, end, code
					new HashMap<String,Map<String,Map<String,Object>>>();
			Map<String, Object> all_res1 = new HashMap<String, Object>();
			Map<String, Object> flags1 = new HashMap<String, Object>();
			JSONOut(String chr, String strand1, int br){
				this.chr = chr; 
				this.strand =strand1;//String.join(";",strand1.split(""));
				all_res1.put("id", expt.get("sampleID"));
				
				
				all_res1.put("id", expt.get("sampleID"));

				
				flags1.put("chr", chr);
				flags1.put("br", br);
				flags1.put("strand", strand);
				String strandc = "L"+strand1.length();
				flags1.put("lens", strandc);
				Set<String> s = Arrays.asList(chr.split(";")).stream().collect(Collectors.toSet());
				Set<String> s1 = s.stream().map(t -> t.substring(0, Math.min(t.length(),3))).collect(Collectors.toSet());
				flags1.put("chroms", s);
//				flags1.put("species", s1);
				String speciesc = "S"+s1.size();
				flags1.put("species", speciesc);
				all_res1.put("flags", flags1);
				all_res1.put("reads", all_res);
				annotate = annotation_mode.contains(strandc) ||
						annotation_mode.contains(speciesc) || 
						annotation_mode.contains("all");
				flags1.put("annot", annotate);
			}
			public  void post() {
				flags1.put("read1", this.first_read);
				flags1.put("nreads", this.reads.size());
			
			//	this.all_res.put("sessionID", sessionID);
				Curl curl = new Curl(all_res1, "addreads");
				Map output = curl.run();
				//executor.execute(curl);
				this.all_res.clear();
				this.reads.clear();
			}
			
			public synchronized void append(String readname,  String endPos, String key,Integer polyA, Integer round){
				 if(reads.size() >0 && (reads.size()+1>=report && !reads.contains(readname))) {
					 this.post();
				 }
				
				if(reads.size()==0) first_read=readname;
				Map<String,Map<String,Object>> all_res_round = all_res.get(round.toString());
				if(all_res_round==null) {
					all_res_round = new HashMap<String,  Map<String, Object>>();
				//	all_res_bin.put(chrom, all_res_chr);
					all_res.put(round.toString(), all_res_round);
				}
				Map<String,Object> all_res_chr_end = all_res_round.get(endPos);
				if(all_res_chr_end==null) {
					all_res_chr_end = new HashMap<String, Object>();
					all_res_round.put(endPos, all_res_chr_end);
				}
				Object vals1 = all_res_chr_end.get(key);
				if(vals1==null) {
					if(!annotate) {
						vals1 = polyA==null ? new int[] {0} : new int [] {0,0,0};
					}else {
						vals1 = new ArrayList<String>();
						//((List<String>)vals1).add(readname);
					}
					all_res_chr_end.put(key, vals1);
				}
				if(vals1 instanceof int[]) {
					int[] vals = (int[])vals1;
					vals[0]+=1;
					if(polyA!=null) {
						if(vals.length==1) {
							vals = new int[] {vals[1],0,0};
						}
						vals[1]+=1;
						vals[2]+=polyA.intValue();
					}
				}else {
					List<String>vals = (List<String>)vals1;
					vals.add(readname);
				}
				reads.add(readname);
			
			}
			
		}
		//int total_num=0;
		public static int report=5000; // how many reads to aggregate before reporting

		
		String first_read = null;
		public Set<String> reads = new TreeSet<String>();
		public synchronized void append(
				StringBuffers sb,
				Integer round){
			Integer polyA = sb.polyA();
			String strand = sb.strand();
			String chrom  = sb.chrom();
			Map<String,JSONOut> all_res_round_strand = all_res.get(strand);
			if(all_res_round_strand==null) {
				all_res_round_strand = new HashMap<String, JSONOut>();
				all_res.put(strand, all_res_round_strand);
			}
			JSONOut all_res_chr = all_res_round_strand.get(chrom);
			if(all_res_chr==null) {
				all_res_chr = new JSONOut(chrom,strand, IdentityProfile1.break_thresh);
				all_res_round_strand.put(chrom, all_res_chr);
			}
			all_res_chr.append(sb.readname, sb.end(round), sb.ref(round), polyA, round);
		}
		
		Map<String, Object> expt;
	
		
		static void printCommand(ProcessBuilder pb) {
			List<String> cmd = pb.command();
			StringBuffer sb  = new StringBuffer();
			for(int i=0; i<cmd.size(); i++){
				sb.append(cmd.get(i)+" ");
			}
			System.err.println(sb.toString());
			
		}
		
		class Curl{
			//String url="http://45.113.234.196:89";
			
				String url1, json_str;
	/* Example			
//				curl -H "Content-Type: application/json" --data '{"sampleID":"results_20240626221127_dorado","flags":{"RNA":true}}' http://45.113.234.196:89/register 
 */
				HttpClient client = HttpClient.newHttpClient();
				Gson gson = new Gson();
				//PrintWriter pw;
				//File endP;
				public Curl(Map map, String endpoint) {
					json_str = gson.toJson(map);
					if(url!=null) {
						url1 = url+"/"+endpoint;
					}
				}
				public Map run() {
					Map output = null;
					if(url==null) {
					//PrintWriter pw = new PrintWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(this.endP, true))));
					if(Outputs.format.equals("json")) {
						Outputs.outputstream.println(json_str);
						Outputs.outputstream.flush();
					}
					}	else {				
					HttpRequest request = HttpRequest.newBuilder()
							.version(Version.HTTP_1_1)
					    .uri(URI.create(url1))
					    .POST(BodyPublishers.ofString(json_str))
					    .setHeader("Content-Type", "application/json")
					    .build();
					
					try {
					HttpResponse<String> response = client.send(request, HttpResponse.BodyHandlers.ofString());
					String str = response.body();
					System.err.println(str);
					output = gson.fromJson(str, Map.class);
					

					} catch (IOException | InterruptedException  e) {
						// TODO Auto-generated catch block
						//e.printStackTrace();
					}
					}
					return output;
				}
				
			}
			
		Number sessionID = null;
		//Executor executor = Executors.newFixedThreadPool(1);
		public void register() {
			if(Outputs.url==null) return; 
			
			Curl curl = new Curl(this.expt, "register");
			Map output = curl.run();
			//all_res2.put("id", this.expt.get("sampleID"));
		
			//all_res2.put("flags", flags2);
		//				executor.execute(curl);
			//			System.err.println(curl.out);
						
		}
		//Map<String, Object> all_res2 = new HashMap<String, Object>();
		
		//Map<String, Object> flags2 = new HashMap<String, Object>();
		
	/*	public Map extract() {
			if(Outputs.url==null) return(null); 
			Curl curl = new Curl(all_res2,"extract");
			Map output = curl.run();
			return output;
		}*/
		
		/** writes the isoform information */
		public synchronized String writeString(CigarCluster cc, int source_index, String chrom, String strand, int round){//, CigarClusters cigarClusters) {
			String id_1 = cc.breaks_hash.secondKey;
			CigarHash2 key = cc.breaks;
			CigarHash2 key2 = cc.cloneBreaks(round);
			List<Integer> startp = cc.start_positions;
				String id_2 = key2.toString(startp, strand);
				String id_ = id_1+"/"+id_2;
				String id = "transcripts/"+id_;
/*				HDFObj obji;
				if(altT!=null){
				if(altT.exists(id)){
				//	m = new HashMap<String, HDFObj>();
					 obji =  new HDFObj(altT.readString(id));
					 obji.expand(new_max_cols);
					 obji.cnts[source_index]+=1;
					// obji.break_sum+=1;
					 for(int j=0; j<key.size(); j++){
						 obji.br[j]+=(float) key.get(j).doubleValue()/CigarCluster.Count.divisor;
					 }
				}else{
					
					obji = new HDFObj( key, source_index, new_max_cols, CigarCluster.Count.divisor); // need to worry about cols? 
				}
				altT.writeString(id,obji.toString());
			}*/
				return id_2;
		}
		
		
		
		private void writeIntMatrix(IHDF5SimpleWriter altT2, String id, int[][] str) {
			if(altT2.exists(id)){
				int[][] matrix = altT2.readIntMatrix(id);
			//	int[][] str1 = merge(matrix, str);
				altT2.writeIntMatrix(id, str);
			}else{
				altT2.writeIntMatrix(id, str);
			}
			
		}

		
		
		public static SAMFileWriter[][] getSamWriter(String chrom,String resdir, String[] in_nmes, boolean presorted, SAMFileHeader header) {
			// TODO Auto-generated method stub
			factB.setUseAsyncIo(true);
		//	factB.
			String[] prefix = "primary:secondary:supplementary:polyA".split(":");
			SAMFileWriter[][] res = new SAMFileWriter[prefix.length][in_nmes.length]; //second row is for supplementary alignments, first for primary alignments
			for(int j=0; j<prefix.length; j++){
				for(int i=0; i<in_nmes.length; i++){
					File f = new File(resdir, in_nmes[i]+"."+chrom+"."+prefix[j]+".bam");
					System.err.println("new sam writer "+f.getAbsolutePath());
					res[j][i]  = factB.makeBAMWriter(header, presorted, f);//factory.newWriter(f);
				}
			}
			return res;
		}

		public static FastqWriter[][] getFqWriter(String chrom,String resdir, String[] in_nmes) {
			// TODO Auto-generated method stub
			String[] prefix = "primary:secondary:supplementary:polyA".split(":");
			FastqWriter[][] res = new FastqWriter[prefix.length][in_nmes.length]; //second row is for supplementary alignments, first for primary alignments
			for(int j=0; j<prefix.length; j++){
				for(int i=0; i<in_nmes.length; i++){
					File f = new File(resdir, in_nmes[i]+"."+chrom+"."+prefix[j]+".fastq");
					System.err.println("new fastq writer "+f.getAbsolutePath());
					res[j][i]  = factory.newWriter(f);
				}
			}
			return res;
		}

		public static void waitOnThreads(ExecutorService executor, int sleep) {
			if(executor==null) return ;
			if(executor instanceof ThreadPoolExecutor){
		    	while(((ThreadPoolExecutor) executor).getActiveCount()>0){
		    		try{
			    	System.err.println("OUTPUTS: awaiting completion "+((ThreadPoolExecutor)executor).getActiveCount());
			    	//Thread.currentThread();
					Thread.sleep(sleep);
		    		}catch(InterruptedException exc){
		    			exc.printStackTrace();
		    		}
		    	}
		    	}
			
		}







//		public  void writeUnmapped(String readName) {
			
	//		unmapped.println(readName);
			
	//	}





 //keep on skippping until you hit the target

		public boolean skip(String readName) {
			if(last_id==null){
			  	return false;
			}else{
				if(readName.equals(last_id)){
					last_id= null;
				}
				return true;
			}
		
		}







		public static PrintStream getOutput(File output_join, boolean append) {
		//	File output_join = resdir==null ? new File(output_join1)  : new File(resdir, output_join1);
			PrintStream p1 = null;
			try {
				p1 = (output_join.getName().endsWith(".gz") ? new PrintStream(new GZIPOutputStream(new FileOutputStream(output_join, append))) : new PrintStream(new FileOutputStream(output_join, append)));
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			return p1;
		}







		public static void makeOutputs(String resDir1, boolean append) throws FileNotFoundException, IOException {
			Outputs.resDir = new File(resDir1);
			resDir.mkdir();
			System.err.println("resDir "+resDir);
			boolean gz = gzip;
			List<String> outps = Arrays.asList("join:overlap:nogap:splice:3:5:all:rem".split(":")); //"all"
			List<File> files = outps.stream().map(t-> new File(resDir1, gz? t+".fa.gz" : t+".fa")).collect(Collectors.toList());
			
			if(resDir.listFiles().length>0 && !append) {
				files.stream().map(f->f.delete());
				//throw new IOException("need to do append if the directory exists, or delete "+resDir1);
			}
		
		
			
			Outputs.ps= files.stream().map(t->getOutput(t, append)).collect(Collectors.toList());
			Outputs.joinOut = ps.get(0);
			Outputs.overlapOut = ps.get(1);
			Outputs.noGap = ps.get(2);
		
			Outputs.spliceOut = ps.get(3);
			Outputs.threeOut = ps.get(4);
			Outputs.fiveOut = ps.get(5);
			Outputs.remOut = ps.get(7);
			Outputs.allOut = ps.get(6);
			//return (output_join1.endsWith(".gz") ? new PrintStream(new GZIPOutputStream(new FileOutputStream(output_join, append))) : new PrintStream(new FileOutputStream(output_join, append)));

			//Outputs.getOutput(resDir, t, append, gz)).collect(Collectors.toList());
			
		}


		static List<PrintStream> ps = new ArrayList<>();



		

		

		

		//public static void shutDownExecutor() {
		//	if(executor!=null) executor.shutdown();
			
	//	}

		

	

		

		
		


		

		
	}