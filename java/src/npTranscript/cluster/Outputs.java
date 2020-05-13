package npTranscript.cluster;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.io.OutputStreamWriter;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.zip.GZIPOutputStream;

import ch.systemsx.cisd.hdf5.HDF5Factory;
import ch.systemsx.cisd.hdf5.IHDF5SimpleWriter;
import htsjdk.samtools.fastq.FastqRecord;
import htsjdk.samtools.fastq.FastqWriter;
import htsjdk.samtools.fastq.FastqWriterFactory;
import japsa.bio.np.ErrorCorrection;
import japsa.seq.Alphabet;
import japsa.seq.Sequence;
import japsa.seq.SequenceOutputStream;
import npTranscript.cluster.CigarCluster.Count;
import npTranscript.run.CompressDir;

public class Outputs{
	static final FastqWriterFactory factory = new FastqWriterFactory();
	
	 class FOutp{
		//String nme;
		//boolean gz;
		SequenceOutputStream os;
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
				os =new SequenceOutputStream(os1);
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
	public static boolean mergeSourceClusters = true;
	public static boolean gzipFasta = false;
	public static boolean keepAlignment = true;
	public static boolean keepinputFasta = true;
	public static boolean writeUnSplicedFastq = false;
	
		public File transcripts_file;
		public File reads_file; 
		private final File outfile, outfile1, outfile2, outfile2_1, outfile4, outfile5, outfile6, outfile7, outfile10;
		//outfile9;
		private final FOutp[] leftover_l, unspliced;//, leftover_r, fusion_l, fusion_r;
	
		final int seqlen;
		 PrintWriter transcriptsP,readClusters;
		 IHDF5SimpleWriter clusterW = null;
		 IHDF5SimpleWriter altT = null;
		File resDir;
		CompressDir[] clusters;
		String[] type_nmes;
		public void close() throws IOException{
			//transcriptsP.close();
			//readClusters.close();
			clusterW.close();
			this.altT.close();
			//this.clusters.close();
			for(int i=0; i<clusters.length; i++){
				//if(so[i]!=null) this.so[i].close();
				if(clusters[i]!=null) this.clusters[i].run();
				
			}
			for(int i=0; i<leftover_l.length; i++){
				if(leftover_l[i]!=null) this.leftover_l[i].close();
				if(unspliced[i]!=null) this.unspliced[i].close();
			}
		}
			
		
		boolean writeDirectToZip = false;
		
		int genome_index=0;
		//final FOutp[] so;
	//	int currentIndex=0;
		
		/*public void updateChromIndex(int currentIndex2) throws IOException{
			this.genome_index = currentIndex2;
			this.newReadCluster(currentIndex2);
		}*/
		/*private void newReadCluster(int genome_index) throws IOException{
			 if(readClusters!=null) readClusters.close();
			
		}*/
		public Outputs(File resDir,  String[] type_nmes, boolean overwrite, int currentIndex, boolean isoforms, boolean cluster_depth, int seqlen) throws IOException{
			this.type_nmes = type_nmes;
			this.genome_index= currentIndex;
			 this.resDir = resDir;
			 this.seqlen = seqlen;
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
		//	String prefix = readsF.getName().split("\\.")[0];
			
			 if(doMSA!=null && mergeSourceClusters){
				 clusters = new CompressDir[] {new CompressDir(new File(resDir,  genome_index+"." +"clusters"))};
		//		 so =  new FOutp[] {new FOutp("consensus")};
			 }else if(doMSA!=null && !mergeSourceClusters){
				 clusters =  new CompressDir[type_nmes.length];
			//	 this.so = new FOutp[type_nmes.length];
				 for(int i=0; i<clusters.length; i++){
					 String nmei =  genome_index+"."+type_nmes[i]+".";
					 clusters[i] = new CompressDir(new File(resDir, nmei+"clusters"));
				//	 so[i] = new FOutp(nmei+"consensus" );
				 }
			 }else{
				 //no msa
				clusters = new CompressDir[1];
			//	so = new FOutp[1];
			 }
			 leftover_l = new FOutp[type_nmes.length];
			 unspliced = new FOutp[type_nmes.length];
		/*	 leftover_r = new FOutp[type_nmes.length];
			 fusion_l = new FOutp[type_nmes.length];
			 fusion_r = new FOutp[type_nmes.length];*/
		//	 String nmei = genome_index+"."
			 for(int i=0; i<leftover_l.length; i++){
				 this.leftover_l[i]=  new FOutp(type_nmes[i]+".leftover", true);
				 if(writeUnSplicedFastq) this.unspliced[i]=  new FOutp(type_nmes[i]+".unspliced", true);
				/* this.leftover_r[i]=  new FOutp(type_nmes[i]+".leftover_r" , true);
				 this.fusion_l[i]=  new FOutp(type_nmes[i]+".fusion_l" , true);
				 this.fusion_r[i]=  new FOutp(type_nmes[i]+".fusion_r" , true);*/
			 }
		//	this.right=  new SequenceOutputStream((new FileOutputStream(new File(resDir,genome_index+".right" ))));
			 reads_file = new File(resDir,genome_index+ ".readToCluster.txt.gz");
			 readClusters = new PrintWriter(
					new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(reads_file))));
			 readClusters.println("readID\tclusterId\tsubID\tsource\tlength\tstart_read\tend_read\t"
			 		+ "type_nme\tchrom\tstartPos\tendPos\tbreakStart\tbreakEnd\terrorRatio\tupstream\tdownstream\tstrand\tbreaks");
		
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

		public PrintWriter getBreakPointPw( String chrom, int i)  throws IOException{
			System.err.println("writing breakpoint files");
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

		
		
		public void writeToCluster(String ID, String subID,  int i, Sequence seq, String str, String name) throws IOException{
			CompressDir cluster = this.getCluster(i);
			seq.setName(name);
			String entryname =  ID;
			if(subID!=null) entryname = entryname+"."+subID;
			if(cluster!=null){
				SequenceOutputStream out   = cluster.getSeqStream(entryname+".fa", true);
				seq.writeFasta(out);
				out.close();
				cluster.closeWriter(out);
				if(str!=null){
					OutputStreamWriter osw = cluster.getWriter(entryname, false, true);
					osw.write(str);osw.write("\n");
					cluster.closeWriter(osw);
				}
			}
		}
		
		public void writeIntMatrix(String id, int[][] matr) {
			this.clusterW.writeIntMatrix(id, matr);
			
		}
		
		private CompressDir getCluster(int source){
			if(clusters==null) return null;
			return mergeSourceClusters? this.clusters[0] : this.clusters[source];
		}
		private int getCount(Count val, int source){
			return mergeSourceClusters ? Arrays.stream(val.count()).sum()  : val.count()[source];
		}

		public void msa(String ID, int source,  CigarCluster cc) throws IOException, InterruptedException{
			// TODO Auto-generated method stub
			if(true) return;
			//NOTE no longer calculating msa
			
			//int source = mergeSourceClusters ? 0 : source1;
			CompressDir cluster = getCluster(source);
			if(cluster==null) return ;
			Iterator<Entry<CigarHash2, Count>> it  = cc.all_breaks.entrySet().iterator();
			for(; it.hasNext();){
					Entry<CigarHash2,Count> ch = it.next();
					CigarHash2 key = ch.getKey();
					Count val = ch.getValue();
					
					String ID1 = ID+"."+val.id();
			File faiFile =new File(cluster.inDir, ID1+".fa");//f[i].getName()+""); 
				if(faiFile.exists() && val.count()[source]>IdentityProfile1.msaDepthThresh){
				File faoFile =new File(cluster.inDir, ID1+".fasta");
	    		ErrorCorrection.runMultipleAlignment(faiFile.getAbsolutePath(),
	    				faoFile.getAbsolutePath());
	    		ArrayList<Sequence> seqList;
	    		if(faoFile.exists()){
	    			seqList = ErrorCorrection.readMultipleAlignment(faoFile.getAbsolutePath());
	    			Sequence consensus = ErrorCorrection.getConsensus(seqList);
					consensus.setDesc("ID="+ID+";subID="+val.id()+";key="+cc.breaks_hash.secondKey+";seqlen="+consensus.length()+";count="+seqList.size()+";type_nme="+cc.getTypeNme(seqlen));
					consensus.setName(ID1);
				//	consensus.print(so[source].os);
					if(!keepAlignment) faoFile.delete();
	    		}
			}
				if(!keepinputFasta) faiFile.delete();
			}
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

		public void writeLeft(Sequence subseq,String baseQ,  boolean negStrand, int source_index)  throws IOException{
			FOutp[] leftover = this.leftover_l;
				//	fusion ? (left ? this.fusion_l : this.fusion_r) : (left ? this.leftover_l : this.leftover_r);
			FastqWriter writer = leftover[source_index].fastq;
			if(negStrand){
				subseq = Alphabet.DNA().complement(subseq);
				baseQ = new StringBuilder(baseQ).reverse().toString();
			}
			 writer.write(new FastqRecord(subseq.getName()+ " "+subseq.getDesc(), new String(subseq.charSequence()), "", baseQ));

		//	subseq.writeFasta(leftover[source_index].os);
		//	leftover[source_index].os.flush();
		}

		public void writeUnspliced(Sequence readseq, String baseQ, boolean negStrand, int source_index) {
			if(unspliced==null || unspliced[source_index]==null) return ;
			FastqWriter writer = unspliced[source_index].fastq;
			if(writer==null) return;
			if(negStrand){
				readseq = Alphabet.DNA().complement(readseq);
				baseQ = new StringBuilder(baseQ).reverse().toString();
			}
			 writer.write(new FastqRecord(readseq.getName()+ " "+readseq.getDesc(), new String(readseq.charSequence()), "", baseQ));

			
		}

		

		
	}