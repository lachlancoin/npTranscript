package npTranscript.cluster;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.SortedMap;
import java.util.TreeMap;

/**
 * @author Lachlan Coin
 *
 */


class SparseArrayVector{
static int zero_=0;
//final boolean end;
	public SparseArrayVector(int len){
		this.len = len;
		this.zero = new int[len];
	//	this.end = end;
	}
	private SortedMap<Integer, int[]> m = new TreeMap<Integer,int[]>();
	
	public void transferFrom(SparseArrayVector mapStart) {
		m.putAll(mapStart.m);
		this.valsum+=mapStart.valsum;
	}
	
	private int valsum=0;
	final int len;
	final int[] zero;
	//double valsumByKey=0;
	/*static Integer zero = 0;
	public void addZero(int pos){
		this.m.put(pos, zero);
	}*/
	
	public int[] get1(Integer position){
		int[] arr = m.get(position);
		if(arr==null){
			arr =new int[len];
			Arrays.fill(arr, 0);
			m.put(position, arr);
		}
		return arr;
	}
	public void addToEntry(Integer position, int src_index, int i) {
		valsum+=i;
	//	if(end && position < 500) {
	//		throw new RuntimeException("!!");
	//	}
		this.get1(position)[src_index]++;
	}
	
	private void addToEntry(Integer position, int[] is) {
		int[] arr = this.get1(position);
		//if(end && position < 500) {
		//	throw new RuntimeException("!!");
		//}
		for(int i=0; i<arr.length; i++){
			arr[i]+=is[i];
		}
		
	}
	
	public String toString(){
		return m.keySet().toString();
	}

	public List<Integer> keys() {
		List<Integer> l= new ArrayList<Integer>(this.m.keySet());
		Collections.sort(l);
		return l;
	}
	/*public List<Integer> keys(double thresh) {
		List<Integer> l= new ArrayList<Integer>();
		for(Iterator<Map.Entry<Integer, Integer>> it = m.entrySet().iterator(); it.hasNext();){
			Map.Entry<Integer, Integer> nxt = it.next();
			if(nxt.getValue()>=thresh){
				l.add(nxt.getKey());
			}
		}
		Collections.sort(l);
		return l;
	}*/
	
	public Iterator<Integer> keyIt(){
		return m.keySet().iterator();
	}

	public int[] get(Integer val) {
		int[] res =  this.m.get(val);
		if(res==null) return zero;
		else return res;
	}
	
	
	

	public void clear() {
		m.clear();
		this.valsum=0;
		
	}
	
	public int[] getDepth(Integer i) {
		int[] val = m.get(i);
		return val==null ? zero: val;
	}
	public int getDepth(int src_index, Integer i) {
		int[] val = m.get(i);
		return val==null ? zero_: val[src_index];
	}
	public Iterator<Integer> tailKeys(Integer st) {
		return m.tailMap(st).keySet().iterator();
	}
	void merge(SparseArrayVector source){
		for(Iterator<Integer> it = source.keyIt();it.hasNext();){
			Integer key = it.next();
			
			
			this.addToEntry(key, source.get(key));
		}
		
			//return this.valsum;
		}

	
	public double valsum() {
	return valsum;
	}
	

	
}