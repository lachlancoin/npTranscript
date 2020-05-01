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

import ch.systemsx.cisd.hdf5.HDF5Factory;
import ch.systemsx.cisd.hdf5.IHDF5SimpleWriter;
import japsa.seq.Sequence;
import japsa.seq.SequenceOutputStream;
import npTranscript.cluster.CigarCluster.Count;
import npTranscript.run.CompressDir;

public class Outputs{
		public File transcripts_file;
		public File reads_file; 
		private final File outfile, outfile1, outfile2, outfile2_1, outfile4, outfile5, outfile6, outfile7, outfile10;
		//outfile9;
	
		
		 private PrintWriter transcriptsP,readClusters;
		 IHDF5SimpleWriter clusterW = null;
		 IHDF5SimpleWriter altT = null;
		File resDir;
		CompressDir[] clusters;
		String[] type_nmes;
		public void close() throws IOException{
			transcriptsP.close();
			readClusters.close();
			clusterW.close();
			this.altT.close();
			//this.clusters.close();
			for(int i=0; i<clusters.length; i++){
				this.clusters[i].runMSA(false);
			}
		}
			
		boolean writeDirectToZip = false;
		
		int genome_index=0;
	//	int currentIndex=0;
		
		/*public void updateChromIndex(int currentIndex2) throws IOException{
			this.genome_index = currentIndex2;
			this.newReadCluster(currentIndex2);
		}*/
		/*private void newReadCluster(int genome_index) throws IOException{
			 if(readClusters!=null) readClusters.close();
			
		}*/
		public Outputs(File resDir,  String[] type_nmes, boolean overwrite, int currentIndex, boolean isoforms, boolean cluster_depth) throws IOException{
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
		//	String prefix = readsF.getName().split("\\.")[0];
			 clusters = new CompressDir[type_nmes.length];
			 for(int i=0; i<type_nmes.length; i++){
				 String nmei = genome_index+"."+type_nmes[i]+".clusters";
				 clusters[i] = new CompressDir(new File(resDir, nmei));
			 }
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

		public void writeToCluster(String entryname, int i, Sequence seq, String str) throws IOException{
			SequenceOutputStream out   = clusters[i].getSeqStream(entryname+".fa", true);
			seq.writeFasta(out);
			out.close();
			clusters[i].closeWriter(out);
			
			OutputStreamWriter osw = clusters[i].getWriter(entryname, false, true);
			osw.write(str);osw.write("\n");
			clusters[i].closeWriter(osw);
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