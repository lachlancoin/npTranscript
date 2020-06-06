package npTranscript.cluster;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
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

		//Name,Type,Minimum,Maximum,Length,Direction,gene
	//	5'UTR,5'UTR,1,265,265,forward,none

		public static String NAstring = "NA";
		public static int tolerance = 10;
		public static int correctionDistLeft = 10000;
		public static int correctionDistRight = 10000;	
	
		
		public void print(PrintWriter pw){
			int len = spliced_count.length;
			pw.print("Gene\tStart\tEnd");
			for(int j=0; j<len; j++){
				pw.print("\tSpliced_"+j);
				pw.print("\tUnspliced_"+j);
			}
			pw.println();
			for(int i=0; i<this.start.size(); i++){
				pw.print(this.genes.get(i)+"\t"+start.get(i)+"\t"+end.get(i));
				for(int j=0; j<len; j++){
					pw.print("\t"+this.spliced_count[j][i]);
					pw.print("\t"+this.unspliced_count[j][i]);
				}
				pw.println();
			}
			pw.flush();
		}
		
		public String nextDownstream(int rightBreak, int chrom_index){
			if(rightBreak<0) return null;
			for(int i=0; i<start.size(); i++){
				if(rightBreak -tolerance <= start.get(i) ){//&& rightBreak < end.get(i)){
					return genes.get(i);
				}
			}
			return "end"+(chrom_index>0 ? "."+chrom_index : "");
		}
		
		public String nextUpstream(int leftBreak, int chrom_index){
			if(leftBreak<0) return  "null";
			for(int i=start.size()-1; i>=0 ;i--){
				if(leftBreak >= start.get(i)){// && leftBreak-tolerance<end.get(i)){
					return genes.get(i);
				}
			}
			return "st"+(chrom_index>0 ? "."+chrom_index : "");
		}
	
    int updateLeft(int en, int refStart){
    	int breaks_out_0= en;
    	for(int i=0; i< this.start.size(); i++){
    		int st_i = start.get(i);
	    	if(en+ tolerance > st_i ){
				if(Math.abs(st_i - en) < correctionDistLeft && st_i >refStart) {
					breaks_out_0 = st_i;
					
				}
			}
    	}
    	return breaks_out_0;
    }
    
    int updateRight( int st, int refEnd){
    	for(int i=0; i< this.start.size(); i++){
	    	int st_i = start.get(i);
	    	if(st - tolerance < st_i) {
				int dist = Math.abs(st_i - st);
				if(dist < correctionDistRight && st_i <refEnd) {
					return st_i;
				}else{
					return st;
				}
			}
    	}
    	return st;
    }
		
	public void getBreaks(int[] breaks_in, int[] breaks_out, int refStart, int refEnd){
		System.arraycopy(breaks_in, 0, breaks_out, 0, 2);
		breaks_out[1] = this.updateRight(breaks_in[1], refEnd);
		breaks_out[0] = this.updateLeft(breaks_in[0], refStart);
	}
	public Annotation(String chrom, int seqlen){
		this.seqlen = seqlen;
		this.chrom = chrom;
	}
	public	Annotation(File f,String chrom, int seqlen, int source_count) throws IOException{
		this(chrom, seqlen);
		BufferedReader br = new BufferedReader(new FileReader(f)) ;
			List<String> header = Arrays.asList(br.readLine().split(","));
			int gene_ind = header.indexOf("gene");
			int start_ind = header.indexOf("Minimum");
			int end_ind = header.indexOf("Maximum");
			String str = "";
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
			mkOrfs(source_count);
			
		}
	
	String chrom;
	 void mkOrfs(int source_count) {
		orfs = new int[start.size()][];
		overlap = new double[start.size()];
		orf_len = new int[start.size()];
		this.unspliced_count = new int[source_count][start.size()];
		this.spliced_count = new int[source_count][start.size()];
		for(int i=0; i<orf_len.length; i++) {
			orf_len[i] = end.get(i) - start.get(i)+1;
			orfs[i] = new int[] {start.get(i),end.get(i)};
		}
		
	}
		private int[][] orfs; 
		public int[][] unspliced_count;
		public int[][] spliced_count;

	public  double[] overlap;
	public  int[] orf_len;
	
		public int seqlen() {
			return seqlen;
		}

		public void adjust3UTR(int seqlen2) {
			if(this.end.get(end.size()-1)>seqlen2 && this.start.get(end.size()-1) < seqlen2){
				end.set(end.size()-1, seqlen2);
			}
			
		}

		public boolean isLeader(int prev_pos) {
			return prev_pos < this.start.get(1)-tolerance;
		}

	
		public void addCount(int i, int source_index, boolean spliced) {
			if(spliced) spliced_count[source_index][i]+=1;
			else unspliced_count[source_index][i]+=1;
			
		}
	}