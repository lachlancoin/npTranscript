package npTranscript.cluster;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
/**
 * @author Lachlan Coin
 *
 */

public class CigarHash2 extends ArrayList<Integer> {
	public CigarHash2(boolean rounded){
		this.rounded = rounded;
	}
	//public static boolean subclusterBasedOnStEnd = false;
	boolean rounded = false;
	public CigarHash2 clone(boolean round, int start, int end){
		if(round && this.rounded) {
			throw new RuntimeException("already rounded");
		}
		CigarHash2 obj =new CigarHash2(round || rounded);
		if(round) obj.addAllR(this,start, end);
		else obj.addAll(this, start);
		return obj;
	}
	public CigarHash2 clone(){
		CigarHash2 obj =new CigarHash2(this.rounded);
		obj.addAll(this);
		return obj;
	}
	
	public static CigarHash2 merge(List<CigarHash2> suppl) {
		boolean rounded = suppl.get(0).rounded;
		CigarHash2 obj =new CigarHash2(rounded);
		for(int j=0; j<suppl.size(); j++){
			if(j>0 && suppl.get(j).rounded!=rounded) throw new RuntimeException("!! inconsistent rounding");
			obj.addAll(suppl.get(j));
		}
		Collections.sort(obj);
		if(obj.size() % 2 !=0) {
			for(int i=0; i<suppl.size();i++){
				if(suppl.get(i).size() %2 !=0) {
					throw new RuntimeException("!!");
				}
			}
		}
		return obj;
	}
	
	private void addAll(CigarHash2 cigarHash2, int start) {
		if(this.rounded!=cigarHash2.rounded) throw new RuntimeException("!!");
		for(int i= start; i<cigarHash2.size(); i++){
			add(cigarHash2.get(i));
		}
	}
	
	public void addAllR(CigarHash2 obj, int start) {
		this.addAllR(obj, start, obj.size());
	}

	public void addAllR(CigarHash2 obj, int start, int end) {
		if(this.size()>0) throw new RuntimeException("not empty");
		if(obj.rounded) throw new RuntimeException("is rounded");
		this.rounded = true;
			for(int i=start; i<end; i++){
				addR(obj.get(i));
			}
		
	}
	
	
	public static int round = 100;
	//static int[] breaks_in = new int[2];
	//static int[] breaks_out = new int[2];
//	public static boolean cluster_by_annotation = true;
	/* (non-Javadoc)
	 * @see npTranscript.cluster.CigHash#toString()
	 */
	@Override
	public String toString(){
		return getString(this);
	}
	public String toString(List<Integer> startp) {
		// TODO Auto-generated method stub
		return getString(this, startp);
	}
	
	
	
	public String toString(int round, int st, int end){
		return getString(this, round, st, end);
	}
	private String getString(CigarHash2 l, List<Integer> startp) {
		StringBuffer sb = new StringBuffer();
		for(int j=0; j<startp.size(); j++){
			int from = startp.get(j);
			int to  = j < startp.size()-1 ? startp.get(j+1) : l.size();
			sb.append(getString(l.subList(from, to)));
			if(j<startp.size()-1) sb.append(";");
		}
		return sb.toString();
	}
	
	public static String getString(List<Integer> l){
		StringBuffer sb = new StringBuffer();
		for(int i=0; i<l.size(); i+=2){
			if(i>0) sb.append(",");
			sb.append(l.get(i)+"_"+l.get(i+1));
		}
		return sb.toString();
	}
	public static String getString(List<Integer> l, int toadd){
		StringBuffer sb = new StringBuffer();
		for(int i=0; i<l.size(); i++){
			if(i>0) sb.append(",");
			sb.append(l.get(i)+toadd);
		}
		return sb.toString();
	}
	public static String getString(List<Integer> l, int mult, int st, int end){
		StringBuffer sb = new StringBuffer();
		for(int i=st; i<end; i++){
			if(i>st) sb.append(",");
			sb.append(l.get(i)*mult);
		}
		return sb.toString();
	}
	
	
	private boolean addR(Integer i){ 
		if(!this.rounded) throw new RuntimeException("!!");
		Integer i1 = TranscriptUtils.round(i, round);
		return (super.add(i1));
	}
	
	public static int  overlap(int start, int end, int st1, int end1){
		return (int) Math.min(Math.min(end-start, end1 - st1),Math.min(end1 - start, end-st1));
	}
	public int overlaps(int st, int length) {
		int end = st + length;
	 for(int i=0; i<this.size(); i+=2){
		 int st1 = this.get(i);
		 int end1 = this.get(i+1);
		 double overlap  = overlap(st,end, st1,end1);
		 if(overlap>=0) {
			 return i;
		 }
	 }
	 return -1;
	}
	public List<Integer> rescale() {
		List<Integer> keyv = this.clone();
		for(int j=0; j<keyv.size(); j++){
			keyv.set(j, keyv.get(j)*round);
		}
		return keyv;
	}
	

	



	

	
	
	/* (non-Javadoc)
	 * @see npTranscript.cluster.CigHash#hashCode()
	 */
	/*@Override
	public int hashCode(){
		return this.toString().hashCode();
	}*/
	/*@Override
	public boolean equals(Object obj){
		return this.hashCode()==obj.hashCode();
		//CigarHash2 ch = (CigarHash2)obj;
		//return this.toString().equals(ch.toString());

		/*if(ch.size()!=this.size()) return false;
		for(int i=0; i<size(); i++){
			if(!get(i).equals(ch.get(i))) return false;
		}
		return true;
	}*/
	/*private void roundBreaks(){
		for(int i =0; i<size(); i++){
   			set(i,  TranscriptUtils.round(get(i), round));
		}
   	}*/
	
	
	/*public void adjustBreaks(Annotation annot) {
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
		
	}*/

	
	}