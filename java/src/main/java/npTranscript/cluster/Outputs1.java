package npTranscript.cluster;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.io.PrintWriter;
import java.util.Date;
import java.util.zip.GZIPOutputStream;

import japsa.seq.SequenceOutputStream;

public class Outputs1 {
	public  PrintWriter transcriptsP, annotP,  featureCP, genesP;//, plusMinus;readClusters,

	public  void printTranscript(String str, String depth_str){
		this.transcriptsP.print(str);
		transcriptsP.println(depth_str);
	}
	public  void printGene(String str, String depth_str){
		this.genesP.print(str);
		genesP.println(depth_str);
	}
	public synchronized void printFC(String str){
		this.featureCP.println(str);
	}
	final String genome_index="";
	public File transcripts_file, genes_file;
	public File feature_counts_file, outfile11;
	public static boolean writeGFF=false;
	 SequenceOutputStream refOut;

	 PrintWriter[] gffW;
	 PrintWriter[] bedW;
	
	 public static int[] gffThresh = null;;
		public static boolean firstIsTranscriptome = false;

		public static int gffThreshTranscriptSum = 0;
		
		public static int maxTranscriptsPerGeneInGFF = Integer.MAX_VALUE;
	 
	public Outputs1(File resDir,  String[] type_nmes) throws IOException{
		 int num_sources = type_nmes.length;
			transcripts_file = new File(resDir,genome_index+ "transcripts.txt.gz");
			genes_file = new File(resDir,genome_index+ "genes.txt.gz");
			 feature_counts_file = new File(resDir,genome_index+ "transcripts.fc.txt.gz");
			 outfile11 = new File(resDir,genome_index+ "annotP.txt.gz");
			 String featureCP_header = "Geneid\tChr\tStart\tEnd\tStrand\tLength";

			
				 featureCP =  new PrintWriter( new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(this.feature_counts_file))));

				 transcriptsP =  new PrintWriter( new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(transcripts_file))));
				 genesP =  new PrintWriter( new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(genes_file))));
					this.annotP =  new PrintWriter( new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(this.outfile11))));
					 outfile11 = new File(resDir, genome_index+"annot.txt.gz");
			
				String transcriptP_header = "ID\tchrom\tstart\tend\tnum_exons\tbreaks\texons\tcountTotal\t"+TranscriptUtils.getString("count", num_sources,true);
				String geneP_header = "ID\tchrom\tstart\tend\tnum_exons\tisoforms\texons\tcountTotal\t"+TranscriptUtils.getString("count", num_sources,true);

				/*if(cluster_depth){
					transcriptP_header = transcriptP_header
							+"\t"+TranscriptUtils.getString("depth", num_sources, true)+"\t"+TranscriptUtils.getString("errors", num_sources, true);
							//+"\t"+TranscriptUtils.getString("error_ratio", num_sources, true);
				}*/
				StringBuffer nme_info = new StringBuffer();
				for(int i=0; i<type_nmes.length; i++) nme_info.append(type_nmes[i]+"\t");
				
				if(transcriptsP!=null) transcriptsP.println("#"+nme_info.toString());
				if(transcriptsP!=null) transcriptsP.println(transcriptP_header);
				if(genesP!=null) genesP.println("#"+nme_info.toString());
				if(genesP!=null) genesP.println(geneP_header);
				if(featureCP!=null) featureCP.println("#npTranscript output");
				if(featureCP!=null) featureCP.println(featureCP_header+"\t"+nme_info.toString());

				 if(writeGFF){
					 bedW = new PrintWriter[type_nmes.length];
					 for(int k=0 ;k<bedW.length; k++){
					
						File  bedoutput = new File(resDir,genome_index+k+".bed.gz");
						 bedW[k] = new PrintWriter(
							new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(bedoutput))));
					 bedW[k].println("track name=\""+type_nmes[k]+"\" description=\""+type_nmes[k]+"\" itemRgb=\"On\" ");
					 }
				 }
				 if(writeGFF){
					 gffW= new PrintWriter[firstIsTranscriptome? 2: 1];
					 for(int i=0; i<gffW.length; i++){
						 gffW[i] = new PrintWriter(
								new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(
										new File(resDir,i+".gff.gz")))));
					 Date d = new Date();
						gffW[i].println("##gff-version 3\n"+
								"#description: \n"+
								"#provider: \n"+
								"#contact: \n"+
								"#format: gff3\n"+
								"#date: "+d.toGMTString()+"\n");
					 }
//						this.refOut =new SequenceOutputStream[Annotation.nmes.length];
						//for(int i=0; i<refOut.length; i++){
							File ref_output = new File(resDir,"ref.fa");
							refOut =new SequenceOutputStream(new FileOutputStream(ref_output));
						//}
				 }
	}
	
	public void close() throws IOException{
		
		if(transcriptsP!=null) transcriptsP.close();
		if(genesP!=null) genesP.close();
		if(featureCP!=null) this.featureCP.close();
		if(annotP!=null) this.annotP.close();
		if(bedW!=null){
			for(int i=0; i<bedW.length; i++) bedW[i].close();
		}
		if(gffW!=null) for(int i=0; i<gffW.length; i++) gffW[i].close();
		if(refOut!=null){
				refOut.close();
		}
	}
	
}
