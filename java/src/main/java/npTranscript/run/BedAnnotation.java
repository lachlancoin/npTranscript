package npTranscript.run;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.SortedMap;
import java.util.TreeMap;
import java.util.stream.Collectors;

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
				if(forward) {
					intervalsF.add(new Interval(st,en,gene)) ;
					}else {
						intervalsR.add(new Interval(st,en,gene));
					}
				
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
	
	@Override
	public int nextUpstream(int leftBreak,  boolean forward, int round){
	//	if(leftBreak<0) return null;
		SortedMap<Integer, Integer> indices = new TreeMap<Integer,Integer>();
		for(int i=start.size()-1; i>=0 ;i--){
			if(!enforceStrand || forward==this.strand.get(i)){
			if(leftBreak + round >= start.get(i) && leftBreak-round<=end.get(i)){// && leftBreak-tolerance<end.get(i)){
				
				indices.put(Math.abs(leftBreak - end.get(i)),i);
			}
			}
		}
		
		
		if(indices.size()>0) return indices.get(indices.firstKey());
		return -1;
	}
	

	@Override
	public int nextDownstream(int rightBreak,  boolean forward, int round){
		//if(rightBreak<0) return null;
		SortedMap<Integer, Integer> indices = new TreeMap<Integer,Integer>();
		for(int i=0; i<end.size(); i++){
			if(!enforceStrand || forward==this.strand.get(i)){
			if(rightBreak -round <= end.get(i) && rightBreak+round>=start.get(i) ){//&& rightBreak < end.get(i)){
				indices.put(Math.abs(rightBreak - start.get(i)),i);
				
			}
			}
		}
		if(indices.size()>0) return indices.get(indices.firstKey());
		return genes.size();
		
		
	}
	

	

}
