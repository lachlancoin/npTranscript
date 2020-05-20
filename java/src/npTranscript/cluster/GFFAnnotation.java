package npTranscript.cluster;

import java.io.FileInputStream;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;
import java.util.Set;

import japsa.seq.JapsaAnnotation;
import japsa.seq.JapsaFeature;

public class GFFAnnotation extends Annotation{

	
	
	/* type is gene
	 * 
	 * */
	
	public static Map<String, JapsaAnnotation> readAnno(String gff, String type,Set<String> chrs) throws IOException{
		
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
	public String nextDownstream(int rightBreak, int chrom_index){
		if(rightBreak<0) return null;
		for(int i=0; i<end.size(); i++){
			if(rightBreak -tolerance <= end.get(i) && rightBreak+tolerance>=start.get(i) ){//&& rightBreak < end.get(i)){
				return genes.get(i);
			}
		}
		return   chrom_index+"."+TranscriptUtils.round(rightBreak,CigarHash2.round);
		
	}
	
	@Override
	public String nextUpstream(int leftBreak, int chrom_index){
		if(leftBreak<0) return null;
		for(int i=start.size()-1; i>=0 ;i--){
			if(leftBreak + tolerance >= start.get(i) && leftBreak-tolerance<=end.get(i)){// && leftBreak-tolerance<end.get(i)){
				return genes.get(i);
			}
		}
		return chrom_index+"."+TranscriptUtils.round(leftBreak,CigarHash2.round);
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
	

	
	public	GFFAnnotation(JapsaAnnotation annot, int seqlen, PrintWriter pw) throws IOException{
		super(annot.getAnnotationID(), seqlen);
		String name = "Name";
		String description="description";
		String id = "ID";
		String genenme = "";
		String descr="";
		String biotype = "biotype";
		String biot;
		pw.println("ID\tName\tdescription\tbiotype");
		System.err.println("num features" +annot.numFeatures());
		for(int i=0; i<annot.numFeatures(); i++){
			JapsaFeature f = annot.getFeature(i);
			this.start.add(f.getStart());
			this.end.add(f.getEnd());
			//System.err.println(f.getDesc());
			String[] desc = f.getDesc().split(";");
			 genenme = "";
			 descr = "";
			 biot="";
			String ID = null;
			for(int j=0; j<desc.length; j++){
				String[] descj = desc[j].split("=");
				if(descj[0].equals(name)){
					genenme = descj[1];
				}else if(descj[0].equals(description)){
					descr = descj[1];
				
			}else if(descj[0].equals(biotype)){
				biot = descj[1];
			}
				if(descj[0].equals(id)){
					ID = descj[1].replace("gene:", ""); //.replace("CDS:", "").
				}
			}
			if(ID==null){
				ID = f.getDesc().substring(0,Math.min(f.getDesc().length(),5));
			}
			String str = ID+"\t"+genenme+"\t"+descr+"\t"+biot;
			pw.println(str);
			//System.err.println(str);
			this.genes.add(ID);
			
		}
		
		super.mkOrfs();
	}
}
