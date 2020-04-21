package npTranscript.cluster;

import java.io.FileInputStream;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import japsa.seq.JapsaAnnotation;
import japsa.seq.JapsaFeature;

public class GFFAnnotation extends Annotation{

	
	public static void main(String[] args){
		try{
		String gff = "Chlorocebus_sabaeus.ChlSab1.1.99.gff3.gz";
		String type = "gene";
	//	List<String> genes = Arrays.asList("Name=ACE2:Name=TMPRSS2".split(":"));
	Map<String, JapsaAnnotation> anno =  readAnno(gff, type);
		int seqlen = 1;
		for(int i=0; i<anno.size(); i++){
			GFFAnnotation g_annot = new GFFAnnotation(anno.get(i), seqlen);
		}
		}catch(Exception exc){
			exc.printStackTrace();
		}
	}
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
	public String nextDownstream(int rightBreak){
		if(rightBreak<0) return null;
		for(int i=0; i<end.size(); i++){
			if(rightBreak -tolerance <= end.get(i) && rightBreak+tolerance>=start.get(i) ){//&& rightBreak < end.get(i)){
				return genes.get(i);
			}
		}
		return null;
	}
	
	@Override
	public String nextUpstream(int leftBreak){
		if(leftBreak<0) return null;
		for(int i=start.size()-1; i>=0 ;i--){
			if(leftBreak + tolerance >= start.get(i) && leftBreak-tolerance<=end.get(i)){// && leftBreak-tolerance<end.get(i)){
				return genes.get(i);
			}
		}
		return null;
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
	
	public	GFFAnnotation(JapsaAnnotation annot, int seqlen) throws IOException{
		super(annot.getAnnotationID(), seqlen);
		String name = "Name";
		String id = "ID";
		for(int i=0; i<annot.numFeatures(); i++){
			JapsaFeature f = annot.getFeature(i);
			this.start.add(f.getStart());
			this.end.add(f.getEnd());
			String[] desc = f.getDesc().split(";");
			String genenme = null;
			String ID = null;
			for(int j=0; j<desc.length; j++){
				String[] descj = desc[j].split("=");
				if(descj[0].equals(name)){
					genenme = descj[1];
				}
				if(descj[0].equals(id)){
					ID = descj[1].replace("gene:", ""); //.replace("CDS:", "").
				}
			}
			if(ID==null){
				ID = f.getDesc().substring(0,Math.min(f.getDesc().length(),5));
			}
			//if(genenme==null){
			//	genenme = ID==null ? f.getDesc().substring(0,Math.min(f.getDesc().length(),5)): ID;
			//}
			this.genes.add(ID);
			
		}
		super.mkOrfs();


	}
}
