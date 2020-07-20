package npTranscript.cluster;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.io.OutputStreamWriter;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.concurrent.ExecutorService;
import java.util.zip.GZIPOutputStream;

import ch.systemsx.cisd.hdf5.HDF5Factory;
import ch.systemsx.cisd.hdf5.IHDF5SimpleWriter;
import htsjdk.samtools.fastq.FastqRecord;
import htsjdk.samtools.fastq.FastqWriter;
import htsjdk.samtools.fastq.FastqWriterFactory;
import japsa.seq.Sequence;
import npTranscript.cluster.CigarCluster.Count;
import npTranscript.run.CompressDir;
import npTranscript.run.SequenceOutputStream1;

public class Outputs{
	
	public static  ExecutorService executor ;
	
	static final FastqWriterFactory factory = new FastqWriterFactory();
	
	 class FOutp{
		//String nme;
		//boolean gz;
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
				f = new File(resDir,genome_index+"."+nme+".fastq");
				fastq = factory.newWriter(f);
			}
			else{
				f = new File(resDir,genome_index+"."+nme+".fasta"+(gz ? ".gz": ""));
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
	}
	// public static boolean MSA_at_cluster = true;
	public static List<String> doMSA = null;//"5_3";
	public static Map<Integer,Integer>msa_sources = new HashMap<Integer, Integer>();
//	public static boolean mergeSourceClusters = true;
	public static boolean gzipFasta = false;
	public static boolean keepAlignment = true;
	public static boolean keepinputFasta = true;
	public static boolean writePolyA = false;
	public static int minClusterEntries = 5;
	public static Collection numBreaks = Arrays.asList(new Integer[0]); // numBreaks for MSA 
	
		public File transcripts_file;
		public File reads_file; 
		private final File outfile, outfile1, outfile2,  outfile4, outfile5, outfile6, outfile7, outfile10, outfile11;
		//outfile9;
		private final FOutp[] leftover_l, polyA;//, leftover_r, fusion_l, fusion_r;
	
	//	final int seqlen;
		 PrintWriter transcriptsP,readClusters, annotP;
		 IHDF5SimpleWriter clusterW = null;
		 IHDF5SimpleWriter altT = null;
		File resDir;
		CompressDir[] clusters;
		String[] type_nmes;
		public void close() throws IOException{
			transcriptsP.close();
			readClusters.close();
			
			if(clusterW!=null) clusterW.close();
			this.altT.close();
			this.annotP.close();
			//this.clusters.close();
			for(int i=0; i<clusters.length; i++){
				//if(so[i]!=null) this.so[i].close();
				if(clusters[i]!=null) this.clusters[i].run(Outputs.minClusterEntries *2);
				
			}
			for(int i=0; i<leftover_l.length; i++){
				if(leftover_l[i]!=null) this.leftover_l[i].close();
				if(polyA[i]!=null) this.polyA[i].close();
			}
		}
			
		
		boolean writeDirectToZip = false;
		
		int genome_index=0;
		
		public Outputs(File resDir,  String[] type_nmes, boolean overwrite, int currentIndex, boolean isoforms, boolean cluster_depth) throws IOException{
			this.type_nmes = type_nmes;
			this.genome_index= currentIndex;
			 this.resDir = resDir;
		//	 this.seqlen = seqlen;
			 int num_sources = type_nmes.length;
			 outfile = new File(resDir,genome_index+ ".txt");
			 outfile1 = new File(resDir, genome_index+ "coref.txt");
			 outfile2 = new File(resDir, genome_index+".clusters.h5");
			
			 outfile4 = new File(resDir,genome_index+ ".exons.txt.gz");
			 outfile5 = new File(resDir,genome_index+ ".clusters.fa.gz");
			 outfile6 = new File(resDir,genome_index+ ".tree.txt.gz");
			 outfile7 = new File(resDir,genome_index+ ".dist.txt.gz");
			 outfile10 = new File(resDir,genome_index+".isoforms.h5");
			 
			 outfile11 = new File(resDir, genome_index+".annot.txt.gz");
		//	String prefix = readsF.getName().split("\\.")[0];
			List<Integer> vals = new ArrayList<Integer>(new HashSet<Integer> (Outputs.msa_sources.values()));
			Collections.sort(vals);
			 //List<String>[] types = new ArrayList<String>[vals.size())]; 
			 if(doMSA!=null && ( Outputs.msa_sources.size()==0)){
				 clusters = new CompressDir[] {new CompressDir(new File(resDir,  genome_index+"." +"clusters"), true)};
		//		 so =  new FOutp[] {new FOutp("consensus")};
			 }else if(doMSA!=null &&  ( Outputs.msa_sources.size()>0)){
				 clusters =  new CompressDir[vals.size()];
			//	 this.so = new FOutp[type_nmes.length];
				
				 for(int i=0; i<clusters.length; i++){
					 
					 String nmei =  genome_index+"."+vals.get(i)+".";
					 clusters[i] = new CompressDir(new File(resDir, nmei+"clusters"), true);
				//	 so[i] = new FOutp(nmei+"consensus" );
				 }
			 }else{
				 //no msa
				clusters = new CompressDir[1];
			//	so = new FOutp[1];
			 }
			 leftover_l = new FOutp[type_nmes.length];
			 polyA = new FOutp[type_nmes.length];
		/*	 leftover_r = new FOutp[type_nmes.length];
			 fusion_l = new FOutp[type_nmes.length];
			 fusion_r = new FOutp[type_nmes.length];*/
		//	 String nmei = genome_index+"."
			 for(int i=0; i<leftover_l.length; i++){
				 this.leftover_l[i]=  new FOutp(type_nmes[i]+".leftover", true);
				 if(writePolyA) this.polyA[i]=  new FOutp(type_nmes[i]+".polyA", true);
				/* this.leftover_r[i]=  new FOutp(type_nmes[i]+".leftover_r" , true);
				 this.fusion_l[i]=  new FOutp(type_nmes[i]+".fusion_l" , true);
				 this.fusion_r[i]=  new FOutp(type_nmes[i]+".fusion_r" , true);*/
			 }
		//	this.right=  new SequenceOutputStream((new FileOutputStream(new File(resDir,genome_index+".right" ))));
			 reads_file = new File(resDir,genome_index+ ".readToCluster.txt.gz");
			 readClusters = new PrintWriter(
					new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(reads_file))));
//			 readID  clusterId       subID   source  length  start_read      end_read   
			 //type_nme        chrom   startPos        endPos  breakStart      breakEnd        errorRatio
			 //upstream        downstream      strand  breaks

			 String header = 
 "readID\tclusterId\tsubID\tsource\tlength\tstart_read\tend_read\ttype_nme\tchrom\tstartPos\tendPos\tstrand\tnum_breaks\tleader_break\terrorRatio\tORFs\tstrand\tbreaks\tspan\tspan_count";
			 readClusters.println(header); //\tbreakStart\tbreakEnd\tbreakStart2\tbreakEnd2\tstrand\tbreaks");
		
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
			 
//						cc.id()+"\t"+chrom+"\t"+cc.start+"\t"+cc.end+"\t"+cc.getTypeNme(seqlen)+"\t"+
//					cc.all_breaks.size()+"\t"+((double)cc.breaks.size()-2.0)/2.0+"\t"+cc.breaks_hash.secondKey+"\t"+
//					cc.totLen+"\t"+cc.readCountSum()+"\t"+read_count
//							+"\t"+cc.getTotDepthSt(true)+"\t"+cc.getTotDepthSt(false)+"\t"+cc.getErrorRatioSt());
			 
			 transcriptsP =  new PrintWriter( new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(transcripts_file))));
				String transcriptP_header = "ID\tchrom\tstart\tend\ttype_nme\tisoforms\tnum_breaks\tleader_break\tORFs\tspan\tspan_length"
					+"\ttotLen\tcountTotal\t"+TranscriptUtils.getString("count", num_sources,true)
				+"\t"+TranscriptUtils.getString("depth", num_sources, true)+"\t"+TranscriptUtils.getString("errors", num_sources, true)
				+"\t"+TranscriptUtils.getString("error_ratio", num_sources, true);
				StringBuffer nme_info = new StringBuffer();
				for(int i=0; i<type_nmes.length; i++) nme_info.append(type_nmes[i]+"\t");
				transcriptsP.println("#"+nme_info.toString());
				transcriptsP.println(transcriptP_header);
				this.annotP =  new PrintWriter( new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(this.outfile11))));
			List<String> str = new ArrayList<String>();
			str.add("subID"); //str.add("ẗype"); 
			for(int j=0; j<num_sources; j++){
				str.add("depth"+j);
			}
			if(isoforms){
			altT = 
					 HDF5Factory.open(outfile10);
			altT.writeStringArray("header", str.toArray(new String[0]));
			}
			str .clear();
			str.add("pos"); //str.add("ẗype"); 
			str.add("base");
			for(int j=0; j<num_sources; j++){
				str.add("depth"+j);
			}
			str.set(0, "pos");
			
			for(int j=0; j<num_sources; j++){
				str.add("errors"+j);
			}
			if(cluster_depth){
			clusterW = 
					 HDF5Factory.open(outfile2);
			clusterW.writeStringArray("header", str.toArray(new String[0]));
			}
		//	clusterW.writeStringArray("header", str.toArray(new String[0]));
		}

		public PrintWriter getBreakPointPw( String chrom, int i, int j)  throws IOException{
			System.err.println("writing breakpoint files");
				File outfile1_ =  new File(resDir,chrom+".breakpoints."+type_nmes[i]+"."+j +".txt.gz") ;
//						new File(resDir,chrom+".breakpoints."+type_nmes[i]+"."+i +".txt.gz");
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

		
		
		public void writeToCluster(String ID, String subID,  int i, Sequence seq, String baseQ,  String str, String name, char strand) throws IOException{
			CompressDir cluster = this.getCluster(i);
			if(strand=='-'){
				seq = TranscriptUtils.revCompl(seq);
				baseQ = new StringBuilder(baseQ).reverse().toString();
			}
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
		
		
		
		public void writeLeft(Sequence subseq,String baseQ,  boolean negStrand, int source_index)  throws IOException{
			FOutp[] leftover = this.leftover_l;
			//	fusion ? (left ? this.fusion_l : this.fusion_r) : (left ? this.leftover_l : this.leftover_r);
			FastqWriter writer = leftover[source_index].fastq;
			writeFastq(writer,subseq, baseQ, negStrand, source_index );
		}
			public void writeFastq(FastqWriter writer, Sequence subseq,String baseQ,  boolean negStrand, int source_index)  throws IOException{
				if(writer==null) return;
				Runnable run = new Runnable(){
					@Override
					public void run() {
						 writer.write(new FastqRecord(subseq.getName()+ " "+subseq.getDesc(), 
								 new String(negStrand ? TranscriptUtils.revCompl(subseq).charSequence(): subseq.charSequence()), "", 
								 negStrand ? new StringBuilder(baseQ).reverse().toString() : baseQ));
						
					}
					
				};
			//	run.run();
				executor.execute(run);
			}

		
		public void writePolyA(Sequence readseq, String nme, String baseQ, boolean negStrand, int source_index)  throws IOException{
			if(polyA==null || polyA[source_index]==null) return ;
			FastqWriter writer = polyA[source_index].fastq;
			writeFastq(writer,readseq, baseQ, negStrand, source_index );
			
		}
		
		public void writeIntMatrix(String id, int[][] matr) {
			this.clusterW.writeIntMatrix(id, matr);
			
		}
		
		private CompressDir getCluster(int source){
			if(clusters==null) return null;
			return this.clusters[this.msa_sources.get(source)];
			//return mergeSourceClusters? this.clusters[0] : this.clusters[source];
		}
		private int getCount(Count val, int source){
			return val.count()[this.msa_sources.get(source)];
		//	return mergeSourceClusters ? Arrays.stream(val.count()).sum()  : val.count()[source];
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