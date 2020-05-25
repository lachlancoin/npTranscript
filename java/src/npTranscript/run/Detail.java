package npTranscript.run;

public class Detail{
	public Integer st,en,span_l, span_r, index;
	Detail(int index){
		this.index = index;
	}
	public String toString(){
		return st+"_"+en+"_"+span_l+"_"+span_r+"_"+index;
	}
	public Integer val(String type){
		if(type.equals("st")){
			return st;
		}else if(type.equals("en")){
			return en;
		}else if(type.equals("index")){
		return index;
	}
		return null;
	}
	
}