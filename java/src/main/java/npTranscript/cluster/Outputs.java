package npTranscript.cluster;

import java.awt.Color;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.OutputStream;
import java.io.OutputStreamWriter;
import java.io.PrintWriter;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.SortedSet;
import java.util.TreeSet;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.ThreadPoolExecutor;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;

import org.apache.commons.math3.linear.SparseRealMatrix;

import ch.systemsx.cisd.hdf5.HDF5Factory;
import ch.systemsx.cisd.hdf5.IHDF5SimpleWriter;
import ch.systemsx.cisd.hdf5.IHDF5Writer;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.fastq.FastqRecord;
import htsjdk.samtools.fastq.FastqWriter;
import htsjdk.samtools.fastq.FastqWriterFactory;
import japsa.seq.Sequence;
import npTranscript.NW.PolyAT;
import npTranscript.cluster.CigarCluster.Count;
import npTranscript.run.Barcodes;
import npTranscript.run.CompressDir;
import npTranscript.run.SequenceOutputStream1;
import npTranscript.run.ViralTranscriptAnalysisCmd2;

public class Outputs{
	public static String readsOutputFile=null;
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
	
		
		public File reads_file; 
		
		private final File  outfile2,  outfile10,  outfile12;
		//private final File[] gff_output;
		//outfile9;
		//private final FOutp[] leftover_l, polyA;//, leftover_r, fusion_l, fusion_r;
	
	//	final int seqlen;
		public  PrintWriter readClusters;//, plusMinus;
		public  PrintWriter unmapped;//, plusMinus;

		 IHDF5SimpleWriter clusterW = null;
		 IHDF5Writer altT = null;
		 IHDF5Writer breakPW = null;
		File resDir;
		CompressDir[] clusters;
		
		
		
		String[] type_nmes;
		public void close() throws IOException{
		//	if(plusMinus!=null) plusMinus.close();
			//IdentityProfileHolder.waitOnThreads(100);
			if(writeCompressDirsExecutor!=null){
				Outputs.waitOnThreads(writeCompressDirsExecutor,100);
				
			}
			/*if(fastQwriter!=null){
				Outputs.waitOnThreads(fastQwriter, 100);
			}*/
			if(h5writer!=null){
				Outputs.waitOnThreads(h5writer, 100);
			}
			readClusters.close();
			unmapped.close();
			if(clusterW!=null){
				writeAllDepths();
				clusterW.close();
			}
			if(altT!=null) this.altT.close();
			//this.clusters.close();
			for(int i=0; i<clusters.length; i++){
				if(clusters[i]!=null) this.clusters[i].run(Outputs.minClusterEntries *2, writeCompressDirsExecutor);
			}
			
			
		

			
		}
		
		
			
		
		boolean writeDirectToZip = false;
		
		String genome_index;
		 int[] col_inds, col_inds_depth;  // this is the col_inds for writing count information in h5 library
		 int new_max_cols, new_max_cols_depth;
		private String last_id;
		public Outputs(File resDir,  String[] type_nmes, boolean isoforms, boolean cluster_depth) throws IOException{
			this.type_nmes = type_nmes;
			genome_index = "0.";
			 this.resDir = resDir;
			 int num_sources = type_nmes.length;
			 outfile2 = new File(resDir, genome_index+"clusters.h5");
			 outfile10 = new File(resDir,genome_index+"isoforms.h5");
			outfile12 = new File(library,genome_index+"breakpoints.h5");
			List<Integer> vals = new ArrayList<Integer>(new HashSet<Integer> (Outputs.msa_sources.values()));
			Collections.sort(vals);
			 if(doMSA!=null && ( Outputs.msa_sources.size()==0)){
				 clusters = new CompressDir[] {new CompressDir(new File(resDir,  genome_index+"clusters"), true)};
			 }else if(doMSA!=null &&  ( Outputs.msa_sources.size()>0)){
				 clusters =  new CompressDir[vals.size()];
				 for(int i=0; i<clusters.length; i++){
					 String nmei =  genome_index+vals.get(i)+".";
					 clusters[i] = new CompressDir(new File(resDir, nmei+"clusters"), true);
				 }
			 }else{
				clusters = new CompressDir[1];
			 }
		
		//	this.right=  new SequenceOutputStream((new FileOutputStream(new File(resDir,genome_index+".right" ))));
			 reads_file = readsOutputFile==null ? new File(resDir,"0.readToCluster.txt.gz") : new File(readsOutputFile);
			 if(!overwrite && reads_file.exists()){
				last_id = transferToNewFile(reads_file);//this makes sure there are no partial lines at the end 
			 }else{
				 last_id = null;
					// TODO Auto-generated method stub
			 }
			 OutputStream os = new FileOutputStream(reads_file, !overwrite);
			// boolean gzip = true;
			 if(reads_file.getName().endsWith(".gz")) os = new GZIPOutputStream(os);
			 readClusters = new PrintWriter(
					new OutputStreamWriter(os));
			 
			 OutputStream os1 = new FileOutputStream(reads_file+".unmapped.gz");
				// boolean gzip = true;
				os1 = new GZIPOutputStream(os1);
				 unmapped = new PrintWriter(
						new OutputStreamWriter(os1));
			 
			
//			 readID  clusterId       subID   source  length  start_read      end_read   
			 //type_nme        chrom   startPos        endPos  breakStart      breakEnd        errorRatio
			 //upstream        downstream      strand  breaks

/*if(IdentityProfile1.trainStrand){
	header = header+ "\tleft10bp\tright10bp\tqual_left10bp\tqual_right10bp\tcorrect_strand";
}*/
			 String str1 = 	 ViralTranscriptAnalysisCmd2.barcodes == null ? "": Barcodes.getHeader();
;
if(last_id==null){
			 readClusters.println(IdentityProfile1.header+"\t"+PolyAT.getHeader()+"\t"+str1); //\tbreakStart\tbreakEnd\tbreakStart2\tbreakEnd2\tstrand\tbreaks");
}
			
				
			List<String> str = new ArrayList<String>();
			//str.add("subID"); //str.add("ẗype"); 
			for(int j=0; j<num_sources; j++){
				str.add(this.type_nmes[j]);
			}
			this.col_inds = new int[this.type_nmes.length];
			if(Outputs.writeH5){
				if(overwrite) outfile10.delete(); // remove existing h5 (temporary)
				altT = HDF5Factory.open(outfile10);
			}
			if(Outputs.calcBreaks && Outputs.writeH5){
				this.breakPW =  HDF5Factory.open(outfile12);

			}
			
			if(altT!=null) {
				writeStringArray(altT, "header", str.toArray(new String[0]),this.col_inds,0);
			}
			this.new_max_cols = Arrays.stream(col_inds) .max() .getAsInt()+1;
				     
				     
					
			
			str .clear();
			str.add("pos"); //str.add("ẗype"); 
		//	str.add("base");
			for(int j=0; j<num_sources; j++){
				str.add(type_nmes[j]);
			}
			if(cluster_depth && writeH5){
				if(overwrite) outfile2.delete();
			clusterW = 	 HDF5Factory.open(outfile2);
			col_inds_depth =new int[num_sources];
			writeStringArray(clusterW, "header", str.toArray(new String[0]),col_inds_depth,1);
			this.new_max_cols_depth = Arrays.stream(col_inds_depth) .max() .getAsInt()+1;
			}
			
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
			
			while((st = br.readLine())!=null){
				String[] str = st.split("\t");
				if(str.length!=head.length){
					break;
				}
				pw.println(st);
				last=str[0];
			}
			pw.close();
			br.close();
			String path = reads_file.getAbsolutePath();
			reads_file.delete();
			Files.move(Paths.get(newF.getAbsolutePath()), Paths.get(path));
			
			return last;
			
		}







		/*public IHDF5SimpleWriter getH5Writer(){
			
			return clusterW;
		}*/

		
		private static void writeStringArray(IHDF5SimpleWriter altT2, String string, String[] array, int[] inds,  int start) {
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
		}
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
			this.readClusters.println(string); this.readClusters.flush();
			
		}
		
		private PrintWriter getBreakPointPw( String chrom, int i, int j)  throws IOException{
			System.err.println("writing breakpoint files");
				File outfile1_ =  new File(resDir,chrom+".breakpoints."+type_nmes[i]+"."+j +".txt.gz") ;
//						new File(resDir,chrom+".breakpoints."+type_nmes[i]+"."+i +".txt.gz");
				PrintWriter pw = new PrintWriter(
					new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(outfile1_, !overwrite))));
			
			return pw;
		}
		
		void printMatrix(SparseRealMatrix cod, SparseVector breakSt2, SparseVector breakEnd2,   int chrom_index, int i, int j, String prefix) throws IOException{
			String id = prefix+"/"+this.type_nmes[i]+"/"+j;
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
			this.breakPW.writeIntMatrix(id, mat);
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
		
		public void writeAllDepths(){
			for(Iterator<String> keys = all_maps.keySet().iterator(); keys.hasNext();){
				String key = keys.next(); 
				Maps map = all_maps.get(key);
				writeDepthH5(key, map);
				keys.remove();	
			}
		}
		
		//	cc.breaks_hash.secondKey;
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
		}
		
	
		
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
		
		public synchronized void writeTotalString(CigarClusters cigarClusters){
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
			
		}
		
		/** writes the isoform information */
		public synchronized String writeString(CigarCluster cc, int source_index, String chrom){//, CigarClusters cigarClusters) {
			String id_1 = cc.breaks_hash.secondKey;
			CigarHash2 key = cc.breaks;
			CigarHash2 key2 = cc.cloneBreaks();
			List<Integer> startp = cc.start_positions;
				String id_2 = key2.toString(startp);
				String id_ = id_1+"/"+id_2;
				String id = "transcripts/"+id_;
				HDFObj obji;
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
			}
				return id_;
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







		public  void writeUnmapped(String readName) {
			unmapped.println(readName);
			
		}







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







		

		

		

		//public static void shutDownExecutor() {
		//	if(executor!=null) executor.shutdown();
			
	//	}

		

	

		

		
		


		

		
	}