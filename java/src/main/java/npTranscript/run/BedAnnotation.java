package npTranscript.run;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Arrays;
import java.util.List;

import npTranscript.cluster.Annotation;

public class BedAnnotation extends Annotation {
static List<String> header = Arrays.asList("chrom:Minimum:Maximum:gene:unknown:Direction".split(":"));
static String split="\t";
	public	BedAnnotation(File f,PrintWriter pw, int source_count) throws IOException{
		//this(chrom, seqlen);
		BufferedReader br = new BufferedReader(new FileReader(f)) ;
			
			int gene_ind = header.indexOf("gene");
			int start_ind = header.indexOf("Minimum");
			int end_ind = header.indexOf("Maximum");
			int direction_ind = header.indexOf("Direction");
			String str = "";
			for(int i=0; (str=br.readLine())!=null; i++){
				String[] line = str.split(split);
				boolean both = false;//line[direction_ind].equals("both");
				boolean forward = line[direction_ind].equals("+");
				String gene = line[gene_ind];	
				
				int st = Integer.parseInt(line[start_ind]);
				int en = Integer.parseInt(line[end_ind]);
				if(pw!=null){
					String str1 = "0\t"+gene+"\t"+gene+"\t"+gene+"\t"+gene+"\t"+gene;
					pw.println(str1);
				}
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
			this.seqlen = this.end.get(end.size()-1);
		}
}
