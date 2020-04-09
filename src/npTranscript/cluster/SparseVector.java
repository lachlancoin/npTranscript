package npTranscript.cluster;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;
import java.util.SortedMap;
import java.util.TreeMap;

/**
 * @author Lachlan Coin
 *
 */


class SparseVector{
	static Integer zero = 0;
	
	private SortedMap<Integer, Integer> m = new TreeMap<Integer, Integer>();
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
	public Iterator<Integer> keyIt(){
		return m.keySet().iterator();
	}

	public Integer get(Integer val) {
		Integer res =  this.m.get(val);
		if(res==null) return 0;
		else return res;
	}
	
	
	public double similarity(SparseVector sv_B) {
		double intersection = 0;
		double Asize = m.size();//valsum;// m.size();
		double BnotA = 0;
		for (Iterator<Integer> it = sv_B.keyIt(); it.hasNext();) {
			if (m.containsKey(it.next())) {
				intersection+=1; //already in union
			} else {
				BnotA +=1;
			}
		}
		double union = Asize + BnotA;
		return ((double) intersection) / union;
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
	
	/*
	public static double similarity(Map<Integer, Integer> map, Map<Integer, Integer> m1) {
		int intersection = 0;
		int union = map.size();
		for (Iterator<Integer> it = m1.keySet().iterator(); it.hasNext();) {
			if (map.containsKey(it.next())) {
				intersection++;

			} else {
				union++;
			}
		}
		return ((double) intersection) / (double) union;
	}

	static int merge(Map<Integer, Integer> target, Map<Integer, Integer> source){
	//	int source_sum = sum(source);
	//	int target_sum = sum(target);
	//	int sum0 = source_sum+target_sum;
		Iterator<Entry<Integer, Integer>> it = source.entrySet().iterator();
	
		while (it.hasNext()) {
			Entry<Integer, Integer> entry = it.next();
			Integer key = entry.getKey();
			int curr = target.containsKey(key) ? target.get(key) : 0;
			int newv = curr + entry.getValue();
			target.put(key, newv);
			
		}
	
		return sum(target);
	}*/
}