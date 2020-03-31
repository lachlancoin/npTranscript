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
		
		int[][] orfs; 
		//Name,Type,Minimum,Maximum,Length,Direction,gene
	//	5'UTR,5'UTR,1,265,265,forward,none

		
		public String nextDownstream(int rightBreak, int tolerance){
			for(int i=0; i<start.size(); i++){
				if(rightBreak -tolerance <= start.get(i) && rightBreak < end.get(i)){
					return genes.get(i);
				}
			}
			return "NA";
		}
		public String nextUpstream(int leftBreak, int tolerance){
			for(int i=start.size()-1; i>=0 ;i--){
				if(leftBreak >= start.get(i) && leftBreak-tolerance<end.get(i)){
					return genes.get(i);
				}
			}
			return "NA";
		}
		
	public	Annotation(File f) throws IOException{
			
			BufferedReader br = new BufferedReader(new FileReader(f)) ;
			List<String> header = Arrays.asList(br.readLine().split(","));
			int gene_ind = header.indexOf("gene");
			int start_ind = header.indexOf("Minimum");
			int end_ind = header.indexOf("Maximum");
			String str = "";
			while((str=br.readLine())!=null){
				String[] line = str.split(",");
				String gene = line[gene_ind];
				if(!gene.equals("none")) {
					genes.add(gene);
					int st = Integer.parseInt(line[start_ind]);
					int en = Integer.parseInt(line[end_ind]);
					start.add(st);
					end.add(en);
				
				}
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