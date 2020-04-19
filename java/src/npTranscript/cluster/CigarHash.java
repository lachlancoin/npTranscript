package npTranscript.cluster;

import java.util.ArrayList;
/**
 * @author Lachlan Coin
 *
 */

public class CigarHash extends ArrayList<Integer> {
	public static int round = 100;
	public static boolean cluster_by_annotation = true;
	/* (non-Javadoc)
	 * @see npTranscript.cluster.CigHash#toString()
	 */
	@Override
	public String toString(){
		StringBuffer sb = new StringBuffer();
		for(int i=0; i<this.size(); i++){
			if(i>0) sb.append(",");
			sb.append(this.get(i));
		}
		return sb.toString();
	}
	
	
	@Override
	public boolean add(Integer i){
		//if(this.contains(i)) throw new RuntimeException("error");
		return (super.add(i));
	}
	/* (non-Javadoc)
	 * @see npTranscript.cluster.CigHash#equals(java.lang.Object)
	 */
	@Override
	public boolean equals(Object o){
		if(cluster_by_annotation) return this.secondKey.equals(((CigarHash)o).secondKey);
		else return super.equals(o);
	}


	
	/* (non-Javadoc)
	 * @see npTranscript.cluster.CigHash#clear()
	 */
	@Override
	 public void clear() {
		secondKey = "";
		super.clear();
	 }
	/* (non-Javadoc)
	 * @see npTranscript.cluster.CigHash#hashCode()
	 */
	@Override
	public int hashCode(){
		if(cluster_by_annotation) return secondKey.hashCode();
		StringBuffer sb = new StringBuffer();
   		for(int i =0; i<size(); i++){
   			sb.append(get(i));
   		}
   		return sb.toString().hashCode();
	}
	private void roundBreaks(){
		for(int i =0; i<size(); i++){
   			set(i,  TranscriptUtils.round(get(i), round));
		}
   	}
	/*@Override
public boolean equals(Object o){
		CigarHash cc = (CigarHash) o;
	 return this.equals(cc.breaks);	
	}
}*/
	static int[] breaks_in = new int[2];
	static int[] breaks_out = new int[2];

	
	
	public void adjustBreaks(Annotation annot) {
		if(annot==null) return ;
		int st = get(0);
		if(get(0)<TranscriptUtils.startThresh) set(0,0);
		int end = get(size()-1);
		for(int i=1; i<size()-1; i+=2){
			breaks_in[0] = get(i);
			breaks_in[1] = get(i+1);
			annot.getBreaks(breaks_in, breaks_out, st, end);
			if(breaks_out[0] > get(i-1)){
				set(i, breaks_out[0]);
			}
			if( breaks_out[1] < get(i+2)){
				set(i+1,breaks_out[1]);
			}
		}
		if(annot.seqlen()-get(size()-1)<TranscriptUtils.endThresh) set(size()-1,annot.seqlen);
		this.roundBreaks();
		
	}

	String secondKey="";
	public void setSecondKey(String string) {
		// TODO Auto-generated method stub
		this.secondKey = string;
	}
	}