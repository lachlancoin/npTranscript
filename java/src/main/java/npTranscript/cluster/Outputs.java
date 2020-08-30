package npTranscript.cluster;

import java.awt.Color;
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
import java.util.Date;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.ThreadPoolExecutor;
import java.util.zip.GZIPOutputStream;

import org.apache.commons.math3.linear.SparseRealMatrix;

import ch.systemsx.cisd.hdf5.HDF5Factory;
import ch.systemsx.cisd.hdf5.IHDF5SimpleWriter;
import htsjdk.samtools.fastq.FastqRecord;
import htsjdk.samtools.fastq.FastqWriter;
import htsjdk.samtools.fastq.FastqWriterFactory;
import japsa.seq.Sequence;
import japsa.seq.SequenceOutputStream;
import npTranscript.cluster.CigarCluster.Count;
import npTranscript.run.CompressDir;
import npTranscript.run.SequenceOutputStream1;

public class Outputs{
	
	//public static  ExecutorService executor ;
	
	public static final FastqWriterFactory factory = new FastqWriterFactory();
	public static int gffThresh = 10;
	
	public static final ExecutorService writeCompressDirsExecutor  = Executors.newSingleThreadExecutor();
	public static final ExecutorService fastQwriter = Executors.newSingleThreadExecutor();
	public static final ExecutorService h5writer = Executors.newSingleThreadExecutor();
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
	}
	// public static boolean MSA_at_cluster = true;
	public static List<String> doMSA = null;//"5_3";
	public static Map<Integer,Integer>msa_sources = new HashMap<Integer, Integer>();
//	public static boolean mergeSourceClusters = true;
	public static boolean gzipFasta = false;
	public static boolean keepAlignment = true;
	public static boolean keepinputFasta = true;
	public static boolean writePolyA = false;
	public static boolean writeBed=false;
	public static boolean writeGFF=true;
	public static int minClusterEntries = 5;
	public static Collection numExonsMSA = Arrays.asList(new Integer[0]); // numBreaks for MSA 
	
		public File transcripts_file;
		public File reads_file; 
		private final File  outfile2,  outfile4, outfile5,  outfile10, outfile11, bedoutput, gff_output;
		//outfile9;
		private final FOutp[] leftover_l, polyA;//, leftover_r, fusion_l, fusion_r;
	
	//	final int seqlen;
		 PrintWriter transcriptsP,readClusters, annotP, bedW, gffW;
		 SequenceOutputStream[] refOut;
		 IHDF5SimpleWriter clusterW = null;
		 IHDF5SimpleWriter altT = null;
		File resDir;
		CompressDir[] clusters;
		
		
		
		String[] type_nmes;
		public void close() throws IOException{
			//IdentityProfileHolder.waitOnThreads(100);
			
			
			transcriptsP.close();
			readClusters.close();
			if(bedW!=null) bedW.close();
			if(gffW!=null) gffW.close();
			if(refOut!=null){
				for(int i=0; i<refOut.length; i++){
					refOut[i].close();
				}
			}
			if(clusterW!=null) clusterW.close();
			this.altT.close();
			this.annotP.close();
			//this.clusters.close();
			for(int i=0; i<clusters.length; i++){
				if(clusters[i]!=null) this.clusters[i].run(Outputs.minClusterEntries *2, writeCompressDirsExecutor);
			}
			
			for(int i=0; i<leftover_l.length; i++){
				if(leftover_l[i]!=null) this.leftover_l[i].close();
				if(polyA[i]!=null) this.polyA[i].close();
			}
			if(writeCompressDirsExecutor!=null){
				Outputs.waitOnThreads(writeCompressDirsExecutor,100);
				
			}
			writeCompressDirsExecutor.shutdown();
			Outputs.fastQwriter.shutdown();
			Outputs.h5writer.shutdown();
		}
			
		
		boolean writeDirectToZip = false;
		
		String genome_index;
		String chrom;
		
		public  void updateChrom(Sequence seq, int currentIndex) {
			String chr = seq.getName();
			// TODO Auto-generated method stub
			this.chrom = chr;
			this.chrom = ((chrom.startsWith("chr") || chrom.startsWith("NC")) ? chrom : "chr"+chrom).split("\\.")[0];
			//this.genome_index = currentIndex;
			if(this.gffW!=null){
				gffW.println("##sequence-region "+chr+" "+0+" "+seq.length());
			}
		}
		
		public Outputs(File resDir,  String[] type_nmes, boolean overwrite,  boolean isoforms, boolean cluster_depth) throws IOException{
			this.type_nmes = type_nmes;
		//	this.genome_index= currentIndex;
			genome_index = "0.";
		
			 this.resDir = resDir;
		//	 this.seqlen = seqlen;
			 int num_sources = type_nmes.length;
		//	 outfile = new File(resDir,genome_index+ ".txt");
		//	 outfile1 = new File(resDir, genome_index+ "coref.txt");
			 outfile2 = new File(resDir, genome_index+"clusters.h5");
			
			 outfile4 = new File(resDir,genome_index+ "exons.txt.gz");
			 outfile5 = new File(resDir,genome_index+ "clusters.fa.gz");
			 outfile10 = new File(resDir,genome_index+"isoforms.h5");
			 bedoutput = new File(resDir,genome_index+"bed.gz");
			gff_output = new File(resDir,genome_index+"gff.gz");
			
			 outfile11 = new File(resDir, genome_index+"annot.txt.gz");
		//	String prefix = readsF.getName().split("\\.")[0];
			List<Integer> vals = new ArrayList<Integer>(new HashSet<Integer> (Outputs.msa_sources.values()));
			Collections.sort(vals);
			 //List<String>[] types = new ArrayList<String>[vals.size())]; 
			 if(doMSA!=null && ( Outputs.msa_sources.size()==0)){
				 clusters = new CompressDir[] {new CompressDir(new File(resDir,  genome_index+"clusters"), true)};
				
		//		 so =  new FOutp[] {new FOutp("consensus")};
			 }else if(doMSA!=null &&  ( Outputs.msa_sources.size()>0)){
				 clusters =  new CompressDir[vals.size()];
			//	 this.so = new FOutp[type_nmes.length];
				 for(int i=0; i<clusters.length; i++){
					 
					 String nmei =  genome_index+vals.get(i)+".";
					// if(TranscriptUtils.coronavirus && vals.size()==type_nmes.length) nmei = genome_index+"."+type_nmes[i]+".";
					 clusters[i] = new CompressDir(new File(resDir, nmei+"clusters"), true);
				//	 so[i] = new FOutp(nmei+"consensus" );
				 }
			 }else{
				clusters = new CompressDir[1];
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
			 reads_file = new File(resDir,genome_index+ "readToCluster.txt.gz");
			 readClusters = new PrintWriter(
					new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(reads_file))));
			 if(writeBed){
				 bedW = new PrintWriter(
						new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(this.bedoutput))));
				 bedW.println("track name=\""+"all"+"\" description=\""+"all"+"\" itemRgb=\"On\" ");
			 }
			 if(writeGFF){
				 gffW= new PrintWriter(
							new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(this.gff_output))));
				 Date d = new Date();
					gffW.println("##gff-version 3\n"+
							"#description: \n"+
							"#provider: \n"+
							"#contact: \n"+
							"#format: gff3\n"+
							"#date: "+d.toGMTString()+"\n");
					this.refOut =new SequenceOutputStream[Annotation.nmes.length];
					for(int i=0; i<refOut.length; i++){
						File ref_output = new File(resDir,Annotation.nmes[i]+".ref.fa");
						refOut[i] =new SequenceOutputStream(new FileOutputStream(ref_output));
					}
			 }
//			 readID  clusterId       subID   source  length  start_read      end_read   
			 //type_nme        chrom   startPos        endPos  breakStart      breakEnd        errorRatio
			 //upstream        downstream      strand  breaks

			 String header = 
 "readID\tclusterId\tsubID\tsource\tlength\tstart_read\tend_read\ttype_nme\tchrom\tstartPos\tendPos\tstrand\tnum_breaks\tleader_break\terrorRatio\tORFs\tstrand\tbreaks\tspan\tspan_count";
			 readClusters.println(header); //\tbreakStart\tbreakEnd\tbreakStart2\tbreakEnd2\tstrand\tbreaks");
		
			 transcripts_file = new File(resDir,genome_index+ "transcripts.txt.gz");
			//	newReadCluster(genome_index);

			 File[] f = new File[] { outfile2,  outfile10, transcripts_file};

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
				String transcriptP_header = "ID\tchrom\tstart\tend\ttype_nme\tisoforms\tnum_exons\tleader_break\tORFs\tspan\tspan_length"
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

		
		
		/*public IHDF5SimpleWriter getH5Writer(){
			
			return clusterW;
		}*/

		
		public synchronized void printTranscript(String str){
			this.transcriptsP.println(str);
		}
		public synchronized void printRead(String string) {
			this.readClusters.println(string);
			
		}
		
		private PrintWriter getBreakPointPw( String chrom, int i, int j)  throws IOException{
			System.err.println("writing breakpoint files");
				File outfile1_ =  new File(resDir,chrom+".breakpoints."+type_nmes[i]+"."+j +".txt.gz") ;
//						new File(resDir,chrom+".breakpoints."+type_nmes[i]+"."+i +".txt.gz");
				PrintWriter pw = new PrintWriter(
					new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(outfile1_))));
			
			return pw;
		}
		
		void printMatrix(SparseRealMatrix cod, SparseVector breakSt2, SparseVector breakEnd2,  int chrom_index, int i, int j) throws IOException{
			
			PrintWriter pw = getBreakPointPw(chrom_index+"",i, j);
				
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
		public synchronized void printBed(List<Integer> breaks, String read_name, char strand, int source, String id, String subid, int num_exons, int span, String span_str){
			if(bedW==null) return;
			int col_id = source % col_len;
		
			int startPos = breaks.get(0)-1;
			int endPos = breaks.get(breaks.size()-1)-1;
			StringBuffer block_sizes = new StringBuffer();
			StringBuffer block_start = new StringBuffer();
			String comma = "";
			
			for(int i=0; i<breaks.size(); i+=2){
				block_start.append(comma+(breaks.get(i)-1-startPos));
				block_sizes.append(comma+(breaks.get(i+1)-breaks.get(i)));
				if(i==0) comma=",";
			}
			bedW.println(chrom+"\t"+startPos+"\t"+endPos+"\t"+read_name+"."+id+"\t"+span+"\t"+strand+"\t"+startPos+"\t"+endPos+"\t"
					+col_str[col_id]+"\t"+num_exons+"\t"+block_sizes.toString()+"\t"+block_start.toString()+"\t"+span_str);

		}

		
		
		public synchronized void writeToCluster(String ID, String subID,  int i, Sequence seq, String baseQ,  String str, String name, char strand) throws IOException{
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
			public synchronized void  writeFastq(FastqWriter writer, Sequence subseq,String baseQ,  boolean negStrand, int source_index)  throws IOException{
				if(writer==null) return;
				Runnable run = new Runnable(){
					@Override
					public void run() {
						 writer.write(new FastqRecord(subseq.getName()+ " "+subseq.getDesc(), 
								 new String(negStrand ? TranscriptUtils.revCompl(subseq).charSequence(): subseq.charSequence()), "", 
								 negStrand ? new StringBuilder(baseQ).reverse().toString() : baseQ));
						
					}
					
				};
				Outputs.fastQwriter.execute(run);
			}

		
		public synchronized void writePolyA(Sequence readseq, String nme, String baseQ, boolean negStrand, int source_index)  throws IOException{
			if(polyA==null || polyA[source_index]==null) return ;
			FastqWriter writer = polyA[source_index].fastq;
			writeFastq(writer,readseq, baseQ, negStrand, source_index );
			
		}
		
		public void writeDepthH5(CigarCluster cc, CigarClusters cigarClusters, Sequence chrom, int chrom_index, int totalDepth) {
			if(clusterW!=null && totalDepth>IdentityProfile1.writeCoverageDepthThresh){
			
				cc.addZeros(cigarClusters.seqlen); 
				int[][] matr =cc.getClusterDepth(cigarClusters.num_sources, chrom);
				clusterW.writeIntMatrix(cc.id(), matr);
				
			}
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

		/** writes the isoform information */
		public synchronized void writeString(String id, Map<CigarHash2, Count> all_breaks, int num_sources) {
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
							altT.writeIntMatrix(id, str);
			
		}
		
		
		
		public void writeIsoforms(CigarCluster cc, CigarClusters cigarClusters, Sequence chrom, int chrom_index
				, int totalDepth){
			String id = cc.id();
			
			int[] rcount=cc.readCount;
			
			
			if(altT!=null){
				if(IdentityProfile1.writeIsoformDepthThresh.length==1){
					if(totalDepth>=IdentityProfile1.writeIsoformDepthThresh[0]){
						writeString(id, cc.all_breaks, cigarClusters.num_sources);
					}
				}else{
					boolean write=false;
					for(int j=0; j<rcount.length; j++){
						if(rcount[j]>=IdentityProfile1.writeIsoformDepthThresh[j]){
							write=true;
						}
					}
					if(write){
						writeString(id, cc.all_breaks, cigarClusters.num_sources);
					}
				}
			}
			
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

		

		//public static void shutDownExecutor() {
		//	if(executor!=null) executor.shutdown();
			
	//	}

		

	

		

		
		


		

		
	}