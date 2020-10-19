package npTranscript.cluster;

import java.util.ArrayList;
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


class SparseVector{
	static Integer zero = 0;
	
	SortedMap<Integer, Integer> m = new TreeMap<Integer, Integer>();
	private int valsum=0;
	//double valsumByKey=0;
	
	public void addZero(int pos){
		this.m.put(pos, zero);
	}
	
	public void addToEntry(Integer position, int i) {
		Integer val = m.get(position);
		valsum+=i;
		//valsumByKey+=position.doubleValue()* (double) i;
		m.put(position, val==null ? i : val + i);
	}
	public String toString(){
		return m.keySet().toString();
	}

	public List<Integer> keys() {
		List<Integer> l= new ArrayList<Integer>(this.m.keySet());
		Collections.sort(l);
		return l;
	}
	public List<Integer> keys(double thresh) {
		List<Integer> l= new ArrayList<Integer>();
		for(Iterator<Map.Entry<Integer, Integer>> it = m.entrySet().iterator(); it.hasNext();){
			Map.Entry<Integer, Integer> nxt = it.next();
			if(nxt.getValue()>=thresh){
				l.add(nxt.getKey());
			}
		}
		Collections.sort(l);
		return l;
	}
	public Iterator<Integer> keyIt(){
		return m.keySet().iterator();
	}

	public Integer get(Integer val) {
		Integer res =  this.m.get(val);
		if(res==null) return 0;
		else return res;
	}
	
	
	

	public void clear() {
		m.clear();
		this.valsum=0;
		
	}
	
	public Integer getDepth(Integer i) {
		Integer val = m.get(i);
		return val==null ? zero: val;
	}
	public Iterator<Integer> tailKeys(Integer st) {
		return m.tailMap(st).keySet().iterator();
	}
	void merge(SparseVector source){
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