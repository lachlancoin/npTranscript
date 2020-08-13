package npTranscript.cluster;

import java.io.FileInputStream;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;
import java.util.SortedSet;

import japsa.seq.JapsaAnnotation;
import japsa.seq.JapsaFeature;

public class GFFAnnotation extends Annotation{

	 List<String> type = new ArrayList<String>();
	
	/* type is gene
	 * 
	 * */
	
	public static Map<String, JapsaAnnotation> readAnno(String gff, String type) throws IOException{
		
		FileInputStream aReader = new FileInputStream(gff);		
		ArrayList<JapsaAnnotation> annos = JapsaAnnotation.readMGFF(aReader,0,0,type);
		aReader.close();
		Map<String, JapsaAnnotation> m = new HashMap<String, JapsaAnnotation>();
		for(int i = annos.size()-1; i>=0; i--){
			if(!annos.get(i).isEmpty()){
				JapsaAnnotation anno = annos.get(i);
				m.put(annos.get(i).getAnnotationID(), anno);
			}
		}
		return m;
	}
	@Override
	public String nextDownstream(int rightBreak, int chrom_index, boolean forward){
		if(rightBreak<0) return null;
		for(int i=0; i<end.size(); i++){
			if(!enforceStrand || forward==this.strand.get(i)){
			if(rightBreak -tolerance <= end.get(i) && rightBreak+tolerance>=start.get(i) ){//&& rightBreak < end.get(i)){
				return genes.get(i);
			}
			}
		}
		return   chrom_index+"."+TranscriptUtils.round(rightBreak,CigarHash2.round);//+(enforceStrand ? (forward ? "+" : "-") : "");
		
	}
	
	@Override
	public String nextUpstream(int leftBreak, int chrom_index, boolean forward){
		if(leftBreak<0) return null;
		for(int i=start.size()-1; i>=0 ;i--){
			if(!enforceStrand || forward==this.strand.get(i)){
			if(leftBreak + tolerance >= start.get(i) && leftBreak-tolerance<=end.get(i)){// && leftBreak-tolerance<end.get(i)){
				return genes.get(i);
			}
			}
		}
		return chrom_index+"."+TranscriptUtils.round(leftBreak,CigarHash2.round);//+(enforceStrand ? (forward ? "+" : "-") : "");
	}
	private int nextUpstreamIndex(int leftBreak, boolean forward){
		if(leftBreak<0) return -1;
		for(int i=start.size()-1; i>=0 ;i--){
			if(!enforceStrand || forward==this.strand.get(i)){
			if(leftBreak + tolerance >= start.get(i) && leftBreak-tolerance<=end.get(i)){// && leftBreak-tolerance<end.get(i)){
				return i;
			}
			}
		}
		return -1;
	}
	
	@Override
	public String getString(Collection<Integer> span,SortedSet<String> obj) {
		if(span.size()==0) return "-";
		
		for(Iterator<Integer> it = span.iterator(); it.hasNext();){
			int ii = it.next();
			String paren =this.parents.get(ii); 
			obj.add(paren==null ? this.genes.get(ii) : paren);
			
		}
		StringBuffer sb = new StringBuffer();
		for(Iterator<String> it = obj.iterator(); it.hasNext();){
			sb.append(it.next());
			if(it.hasNext())sb.append(";");
		}
		return sb.toString();
	}
	
	private int nextDownstreamIndex(int rightBreak, boolean forward){
		if(rightBreak<0) return -1;
		for(int i=0; i<end.size(); i++){
			if(!enforceStrand || forward==this.strand.get(i)){
			if(rightBreak -tolerance <= end.get(i) && rightBreak+tolerance>=start.get(i) ){//&& rightBreak < end.get(i)){
				return i;
			}
			}
		}
		return   -1;
		
	}
	
	@Override
	 int updateRight( int st, int refEnd){
	    	for(int i=0; i< this.start.size(); i++){
		    	int end_i = end.get(i);
		    	if(st - tolerance < end_i) {
					int dist = Math.abs(end_i - st);
					if(dist < correctionDistRight && end_i <refEnd) {
						return end_i;
					}else{
						return st;
					}
				}
	    	}
	    	return st;
	    }
	static int tolerance =10;
	static int tolerance1 = 5;
	static Comparator<Entry<String, Set<Integer>>> c = new Comparator<Entry<String, Set<Integer>>>(){

		@Override
		public int compare(Entry<String, Set<Integer>> o1, Entry<String, Set<Integer>> o2) {
			return Integer.compare(o2.getValue().size(),o1.getValue().size());
		}
		
	};
	Map<String, Set<Integer>> m= new HashMap<String, Set<Integer>>();
	Set<Integer> breaks_in_exons = new HashSet<Integer>();
	@Override
	public synchronized String  getSpan(List<Integer> breaks,boolean forward,  Collection<Integer> l, SortedSet<String> geneNames){
		int istart =-1;
		int iend = -1;
		int start_ = breaks.get(0);
		int end_ =breaks.get(breaks.size()-1);
		for(int i=0; i<this.genes.size(); i++){
			if(!enforceStrand || forward==this.strand.get(i)){
				if(this.end.get(i)> start_-tolerance1){
					istart = i;
					break;
				}
			}
		}
		for(int i=this.genes.size()-1;i>=0; i--){
			if(!enforceStrand || forward==this.strand.get(i)){
				if(this.start.get(i)< end_+tolerance1){
					iend= i;
					break;
				}
			}
		}
		if(istart>=0 && iend>=0 && istart <= iend){
			//
			
			for(int i=istart; i<=iend ; i++){
				
				if(span_only.size()==0|| span_only.contains(this.type.get(i))){
					for(int j=0; j<breaks.size(); j++){
						int pos = breaks.get(j);
						if(start.get(i)<=pos+tolerance1 && end.get(i)>=pos - tolerance1){
							
							String paren = this.parents.get(i);
							if(paren==null) paren = genes.get(i);
							Set<Integer>s = m.get(paren);
							if(s==null) m.put(paren, s = new HashSet<Integer>());
							s.add(j);
							l.add(i);
							breaks_in_exons.add(j);
						}	
					}
				}
			}
			
			// this attempting to find the smallest set of genes to explain the breakpoints.  Its not exhaustive search but starts with biggest set first
			if(m.size()>0){
				List<Entry<String, Set<Integer>>> vals =  new ArrayList<Entry<String, Set<Integer>>>(m.entrySet());
				Collections.sort(vals, c);
				Set<Integer> set = vals.get(0).getValue();
				int i=0;
			//	l.add(0);
				geneNames.add(vals.get(0).getKey());
				while(set.size()<breaks_in_exons.size()
						){
					i= i+1;
					set.addAll(vals.get(i).getValue());
					
					geneNames.add(vals.get(i).getKey());
				}
				m.clear();
				breaks_in_exons.clear();
			}
			
			
			StringBuffer sb = new StringBuffer();
			Iterator<String> it = geneNames.iterator();
			while(it.hasNext()){
				sb.append(it.next());
				if(it.hasNext()) sb.append(";");
			}
			return sb.toString();
		}
		return ".";
		
	}
	 public static List<String> span_only = Arrays.asList("protein_coding".split(":"));
	
	@Override
	public String getTypeNme(int start_, int end_, boolean forward) {
		int iend = nextDownstreamIndex(end_, forward);
		int istart = nextUpstreamIndex(start_, forward);
		String left = istart<0 ? "NA" : (start_<= this.start.get(istart)+ TranscriptUtils.startThresh ? "5" : "no5");
		String right = iend<0 ? "NA" : (end_>= this.end.get(iend)- TranscriptUtils.endThresh ? "3" : "no3");
		return left+"_"+right;
		//return null;
	}
	public static void setGFFFeatureNames(String[] str){
		
		name = str[0];
		description = str[1];
		id = str[2];
		biotype = str[3];
		parent = str[4];
	}
	//ID:description:ID:gene_type:gene_name
	
	//ID=exon:ENST00000456328.2:1;Parent=ENST00000456328.2;gene_id=ENSG00000223972.5;transcript_id=ENST00000456328.2;gene_type=transcribed_unprocessed_pseudogene;gene_name=DDX
	
//   ID=exon:ENST00000495576.1:1;Parent=ENST00000495576.1;gene_id=ENSG00000239945.1;transcript_id=ENST00000495

	
	//--GFF_features=exon_id:description:ID:gene_type:gene_name 
//	gene_name:description:ID:gene_type
	
	List<String> parents = new ArrayList<String>();
	public static String name = "Name";
	public static String description="description";
	public static String id = "ID";
	public static String biotype = "biotype";
	public static String parent = "Parent";
	static String emptyString = "";
	public	GFFAnnotation(String chr, JapsaAnnotation annot, int seqlen, PrintWriter pw, int source_count) throws IOException{
		super(annot.getAnnotationID(), seqlen);
	
		String genenme = "";
		String descr="";
		
		//biot;
	//	pw.println("ID\tName\tdescription\tbiotype");
		System.err.println("num features" +annot.numFeatures());
		Map<String, String> biotypes = new HashMap<String,String>();
		for(int i=0; i<annot.numFeatures(); i++){
			JapsaFeature f = annot.getFeature(i);
			String type = f.getType();
			this.start.add(f.getStart());
			this.end.add(f.getEnd());
			this.strand.add(f.getStrand()=='+');
			//System.err.println(f.getDesc());
			String[] desc = f.getDesc().split(";");
			 genenme = "";
			 descr = "";
			 String  biot=null;
			 String paren=null;
			String ID = null;
			for(int j=0; j<desc.length; j++){
				String[] descj = desc[j].split("=");
		//		descj[1] = descj[1].replace("gene:", "").replace("exon:", "");
				String[] descjv = descj[1].split(":");
				String val = descjv[descjv.length-1];
				if(descj[0].equals(name)){
					genenme = val;
				}else if(descj[0].equals(description)){
					descr = val;
				
			}else if(descj[0].equals(biotype)){
				biot = val;
			}else if(descj[0].equals(parent)){
				paren = val;
				paren = paren.replace("transcript:", "");
			}
				
				if(descj[0].equals(id)){
					ID = val; //.replace("CDS:", "").
				}
			}
			if(ID==null){
				ID = genenme;//if.getDesc().substring(0,Math.min(f.getDesc().length(),5));
			}
			String str = chr+"\t"+ID+"\t"+genenme+"\t"+descr+"\t"+biot;
			pw.println(str);
			//System.err.println(str);
			this.genes.add(genenme.length()>0 ? genenme : ID);
			if(biot!=null) type  = biot;
			biotypes.putIfAbsent(type, type);
			this.type.add(biotypes.get(type));
			this.parents.add(paren);
		}
		
		super.mkOrfs(source_count);
	}
	public GFFAnnotation(String chrom, int seqlen) {
		// TODO Auto-generated constructor stub
		super(chrom,seqlen);
	}
}
