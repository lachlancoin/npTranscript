package npTranscript.cluster;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Iterator;
import java.util.List;
import java.util.SortedSet;
/**
 * @author Lachlan Coin
 *
 */

public class Annotation{
		List<String> genes= new ArrayList<String>();
		List<Integer> start = new ArrayList<Integer>();
		 List<Integer> end = new ArrayList<Integer>();
		 List<Boolean> strand = new ArrayList<Boolean>();
		 
		
		 
	//	List<Integer> breakSt =  new ArrayList<Integer>();
		//List<Integer> breakEnd =  new ArrayList<Integer>();
		 final int seqlen;

		//Name,Type,Minimum,Maximum,Length,Direction,gene
	//	5'UTR,5'UTR,1,265,265,forward,none

		public static String NAstring = "NA";
		public static int tolerance = 10;
		public static int correctionDistLeft = 10000;
		public static int correctionDistRight = 10000;	
		public static boolean enforceStrand = true;
		
		public void print(PrintWriter pw){
			
			int len = spliced_count.length;
			pw.print("Gene\tStart\tEnd");
			for(int j=0; j<len; j++){
				pw.print("\tSpliced_5_"+j);
				pw.print("\tUnspliced_5_"+j);
			}
			pw.println();
			for(int i=0; i<start.size(); i++){
				pw.print(genes.get(i)+"\t"+start.get(i)+"\t"+end.get(i));
				for(int j=0; j<len; j++){
					pw.print("\t"+spliced_count[j][i]);
					pw.print("\t"+unspliced_count[j][i]);
				}
				pw.println();
			}
			pw.flush();
				}
				
			
		
		public String nextDownstream(int rightBreak, int chrom_index, boolean forward){
			if(rightBreak<0) return null;
			for(int i=0; i<start.size(); i++){
				if(!enforceStrand || forward==this.strand.get(i)){
					if(rightBreak -tolerance <= start.get(i) ){//&& rightBreak < end.get(i)){
						return genes.get(i);
					}
				}
			}
			return "end"+(chrom_index>0 ? "."+chrom_index : "");//+ (enforceStrand ? (forward ? "+" : "-") : "");
		}
		/*public String convert(List<Integer> l, int chrom_index){
			StringBuffer sb = new StringBuffer();
			for(int i=0; i<l.size(); i++){
				if(l.get(i)==null) sb.append("null");
				else if(l.get(i)<0) sb.append("end"+(chrom_index>0 ? "."+chrom_index : ""));//+ (enforceStrand ? (forward ? "+" : "-") : "");)
				else if(l.get(i)>=start.size()) sb.append("start"+(chrom_index>0 ? "."+chrom_index : ""));//+ (enforceStrand ? (forward ? "+" : "-") : "");)
				else sb.append(genes.get(i));
				if(i%2 ==0){
					sb.append("_") ;
				}else sb.append(","); 
			}
			return sb.toString();
		}*/
		
		public String nextUpstream(int leftBreak, int chrom_index, boolean forward){
			if(leftBreak<0) return  "null";
			for(int i=start.size()-1; i>=0 ;i--){
				if(!enforceStrand || forward==this.strand.get(i)){
				if(leftBreak+tolerance >= start.get(i)){// && leftBreak-tolerance<end.get(i)){
					return genes.get(i);
				}
				}
			}
			return "st"+(chrom_index>0 ? "."+chrom_index : "");//+(enforceStrand ? (forward ? "+" : "-") : "");
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
	
	public	Annotation(File f,String chrom, int seqlen,  PrintWriter pw,int source_count) throws IOException{
		this(chrom, seqlen);
		BufferedReader br = new BufferedReader(new FileReader(f)) ;
			List<String> header = Arrays.asList(br.readLine().split(","));
			int gene_ind = header.indexOf("gene");
			int start_ind = header.indexOf("Minimum");
			int end_ind = header.indexOf("Maximum");
			int direction_ind = header.indexOf("Direction");
			String str = "";
			for(int i=0; (str=br.readLine())!=null; i++){
				String[] line = str.split(",");
				boolean both = line[direction_ind].equals("both");
				boolean forward = line[direction_ind].equals("forward");
				String gene = line[gene_ind];	
				
				int st = Integer.parseInt(line[start_ind]);
				int en = Integer.parseInt(line[end_ind]);
				String str1 = "0\t"+gene+"\t"+gene+"\t"+gene+"\t"+gene+"\t"+gene;
				pw.println(str1);
				genes.add(gene);
				start.add(st);
				end.add(en);
				if(both && enforceStrand){
					strand.add(true); strand.add(false);
					genes.add(gene);
					start.add(st);
					end.add(en);
				}else{
					strand.add(forward);
				}
			
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

		public synchronized void adjust3UTR(int seqlen2) {
			if(end.size()==0) return ;
			if(this.end.get(end.size()-1)>seqlen2 && this.start.get(end.size()-1) < seqlen2){
				System.err.println("adjusting 3UTR end "+seqlen2);
				end.set(end.size()-1, seqlen2);
			}
			
		}

		public boolean isLeader(int prev_pos) {
			return prev_pos < this.start.get(1)-tolerance;
		}

	
		public void addCount(int i, int source_index, boolean spliced) {
			if(spliced) {
				spliced_count[source_index][i]+=1;
				//if(i==2){
				//	System.err.println("hh");
			//	}
			}
			else unspliced_count[source_index][i]+=1;
			
		}
		public String getTypeNme(int start, int end, boolean forward) {
			return nmes[getTypeInd(start,end, forward)];
			//return null;
		}
		public int getTypeInd(int start, int end, boolean forward) {
			if(start <=TranscriptUtils.startThresh) return end >= seqlen -TranscriptUtils.endThresh ? 0:1;
			else return end >= seqlen -TranscriptUtils.endThresh ? 2:3;
			//return null;
		}

		public static String[] nmes = new String[] {"5_3", "5_no3", "no5_3", "no5_no3"};
		public  char getStrand(Iterator<Integer> l){
			if(!l.hasNext()) return 'N';
			Boolean strandi = this.strand.get(l.next());
			while(l.hasNext()){
				if(!strand.get(l.next()).equals(strandi)){
					return 'N';
				}
			}
			return strandi ? '+' : '-';
			
		}

		public String getSpan(List<Integer>breaks, boolean forward, Collection<Integer> span, SortedSet<String> parents) {
			// TODO Auto-generated method stub
			return ".";
		}

		public String getString(Collection<Integer> span, SortedSet<String> geneNmes) {
			if(span.size()==0) return "-";
			StringBuffer sb = new StringBuffer();
			for(Iterator<Integer> it = span.iterator(); it.hasNext();){
				String gene = this.genes.get(it.next());
				geneNmes.add(gene);
				sb.append(gene);
				if(it.hasNext())sb.append(";");
			}
			return sb.toString();
		}
	}