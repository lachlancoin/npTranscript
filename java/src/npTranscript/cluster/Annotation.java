package npTranscript.cluster;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import npTranscript.cluster.TranscriptUtils;
/**
 * @author Lachlan Coin
 *
 */

public class Annotation{
		List<String> genes= new ArrayList<String>();
		List<Integer> start = new ArrayList<Integer>();
		List<Integer> end = new ArrayList<Integer>();
	//	List<Integer> breakSt =  new ArrayList<Integer>();
		//List<Integer> breakEnd =  new ArrayList<Integer>();
final int seqlen;
		int[][] orfs; 
		//Name,Type,Minimum,Maximum,Length,Direction,gene
	//	5'UTR,5'UTR,1,265,265,forward,none

		static String NAstring = "NA";
		public String nextDownstream(int rightBreak){
			if(rightBreak<0) return NAstring;
			for(int i=0; i<start.size(); i++){
				if(rightBreak -tolerance <= start.get(i) ){//&& rightBreak < end.get(i)){
					return genes.get(i);
				}
			}
			return NAstring;
		}
		public String nextUpstream(int leftBreak){
			if(leftBreak<0) return NAstring;
			for(int i=start.size()-1; i>=0 ;i--){
				if(leftBreak >= start.get(i)){// && leftBreak-tolerance<end.get(i)){
					return genes.get(i);
				}
			}
			return NAstring;
		}
	static int tolerance = 10;
	 static int correctionDistLeft = 10000;
	static int correctionDistRight = 10000;

	public void getBreaks(int[] breaks_in, int[] breaks_out, int refStart, int refEnd){
		int st = breaks_in[1]; //start after break
		int en = breaks_in[0];
		System.arraycopy(breaks_in, 0, breaks_out, 0, 2);
		for(int i=0; i< this.start.size(); i++){
			int st_i = start.get(i);
			if(st - tolerance < st_i) {
				int dist = Math.abs(st_i - st);
				if(dist < correctionDistRight && st_i <refEnd) {
					breaks_out[1] = st_i;
				}
				break;
			}
			if(en+ tolerance > st_i ){
				if(Math.abs(st_i - en) < correctionDistLeft && st_i >refStart) {
					breaks_out[0] = st_i;
				}
			}
		}
		
	}
	public	Annotation(File f, int seqlen) throws IOException{
		this.seqlen = seqlen;
			BufferedReader br = new BufferedReader(new FileReader(f)) ;
			List<String> header = Arrays.asList(br.readLine().split(","));
			int gene_ind = header.indexOf("gene");
			int start_ind = header.indexOf("Minimum");
			int end_ind = header.indexOf("Maximum");
			String str = "";
			//int break_start = -1;
			for(int i=0; (str=br.readLine())!=null; i++){
				String[] line = str.split(",");
				String gene = line[gene_ind];	
					genes.add(gene);
					int st = Integer.parseInt(line[start_ind]);
					int en = Integer.parseInt(line[end_ind]);
					start.add(st);
					end.add(en);
				//	if(i>0){
				//		this.breakSt.add(break_start);
				//		this.breakEnd.add(st);
				//	}
					//break_start = st;
			
			}
			br.close();
			orfs = new int[start.size()][];
			overlap = new double[start.size()];
			orf_len = new int[start.size()];
			
			for(int i=0; i<orf_len.length; i++) {
				orf_len[i] = end.get(i) - start.get(i)+1;
				orfs[i] = new int[] {start.get(i),end.get(i)};
			}
			
		}
	/*public String getInfo(int i) {
			if(i<0 || overlap[i] <0) return "NA";
			else return this.genes.get(i)+","+this.start_offset[i]+","+this.end_offset[i];
		}*/
	public final double[] overlap;
	public final int[] orf_len;
	
	
	
	//public final int[] start_offset; //relative start position of ORF
//	public final int[] end_offset; //relative end position of ORF
		public String calcORFOverlap(int start, int end, int[] first_last, int transcript_len) {
			//could consider correcting here to keep in-fram
			int first = -1;
			int last = -1;
			StringBuffer sb = new StringBuffer();
			for(int i=0 ; i<orfs.length; i++) {
				int[] st_en = orfs[i];
			
				overlap[i] =  (double) TranscriptUtils.overlap(start, end, st_en[0], st_en[1])/ (double)(st_en[1] - st_en[0] +1);
				if(overlap[i]>0) {
					if(first<0) first = i;
					last = i;
					sb.append(";");
					sb.append(this.genes.get(i));
					sb.append(",");
					sb.append(String.format( "%5.3g",overlap[i]).trim());
					sb.append(",");
					sb.append(st_en[0] - start + transcript_len); // how far into read ORF starts;
					sb.append(",");
					sb.append(st_en[1] - start+transcript_len); // how far into read ORF ends
				}
			}
			first_last[0] =first;
			first_last[1] = last;
			return sb.toString();
		}
	}