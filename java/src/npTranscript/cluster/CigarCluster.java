package npTranscript.cluster;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;

import japsa.seq.Sequence;


/**
 * @author Lachlan Coin
 *
 */

public class CigarCluster  {
		//final int index;
		
		
		static int round2 = 100;
		public static boolean recordDepthByPosition = false;
		int breakSt = -1;
		int breakEnd = -1;
		int breakSt2 = -1;
		int breakEnd2 = -1;
		 Boolean forward = null;
		
		private final String id;
		
		public String id(){
			return id;
		}
		/*public String id(int chrom_index){
			return chrom_index+"."+id;
		}*/

		 int start=0;
		 int end=0;
		int prev_position =0;
		int break_point_cluster = -1;		
		boolean first = true;
		
	
		public void addReadCount(int source_index) {
			readCount[source_index]++;
			this.readCountSum++;
			
		}

	public static class Count{
		public Count(int[] count, int id){
			this.count = count;
			this.id = id;
		}
		public Count(int num_sources, int src_index, int id) {
			this.count = new int[num_sources];
			count[src_index]=1;
			this.id = id;
			// TODO Auto-generated constructor stub
		}
		private int id;
		private int[] count;
		public void increment(int source_index) {
			this.count[source_index] = this.count[source_index]+1;
			
		}
		public int[] count(){
			return count;
		}
		public int id(){
			return id;
		}
		
	}
		
   	Map<CigarHash2, Count> all_breaks  ;
   	CigarHash2 breaks = new CigarHash2();
   	CigarHash breaks_hash = new CigarHash();
   	
   	

	
		
		/*public void setBreaks(CigarHash2 breaks){
			this.breaks.addAll(breaks);
		}*/

		public CigarCluster(String id,  int num_sources){
			this.id = id;
			this.readCount = new int[num_sources];
			if(recordDepthByPosition){
				map = new SparseVector();
				this.maps = new SparseVector[num_sources];
				this.errors= new SparseVector[num_sources];
				for(int i=0; i<maps.length; i++){
					maps[i] = new SparseVector();
					errors[i] = new SparseVector();
					
				}
			}else{
				maps = null;
				errors = null;
				map = null;
			}
		}
		
		public CigarCluster(String id,  int num_sources, CigarCluster c1, int source_index) throws NumberFormatException{
			this(id, num_sources);
			this.span.addAll(c1.span);
			this.forward = c1.forward;
			this.breakSt = c1.breakSt;
			this.breakEnd = c1.breakEnd;
			this.breakSt2 = c1.breakSt2;
			this.breakEnd2 = c1.breakEnd2;
			this.breaks.addAllR(c1.breaks);
			all_breaks = new HashMap<CigarHash2, Count>();
			//if(all_breaks.size()>0) throw new RuntimeException("should be zero");
			
			this.all_breaks.put(breaks,new Count(num_sources, source_index,  0));
			this.breaks_hash.secondKey = c1.breaks_hash.secondKey;
			addReadCount(source_index);
			start = c1.start;
			end = c1.end;
			if(recordDepthByPosition){
				map.merge( c1.map);
				for(int i=0; i<maps.length; i++){
					maps[i].merge( c1.maps[i]);
					errors[i].merge(c1.errors[i]);
					if(errors[i].valsum()<0) throw new NumberFormatException("shoud not be negative");
					if(errors[i].valsum()>maps[i].valsum()) throw new NumberFormatException("shoud not be greater");
				}
			}
		}

		private final SparseVector map;// = new SparseVector(); //coverage at high res
		final private SparseVector[] maps, errors;

		public void clear(int source_index) {
			
			
			span.clear();
			this.forward = null;
			this.breakSt=-1;
			this.breakEnd = -1;
			this.breakSt2=-1;
			this.breakEnd2 = -1;
			if(this.recordDepthByPosition){
				map.clear();
				for(int i=0; i<maps.length; i++){
					maps[i].clear();
					errors[i].clear();
				}
			}
			Arrays.fill(readCount, 0);
			readCount[source_index]=1;
			readCountSum=1;
			
			this.prev_position = 0;
			this.break_point_cluster = -1;
			first = true;
			this.breaks.clear();			
			if(all_breaks!=null)this.all_breaks.clear();
			this.breaks_hash.clear();
		}
	

		public void add(int pos, int src_index, boolean match) {
			if(match){
			//System.err.println(prev_position + " "+pos);
			if(first){
				first = false; 
				breaks.add(pos);
			} else{	
				if(pos-prev_position>TranscriptUtils.break_thresh){
					this.breaks.add(prev_position);
					this.breaks.add(pos);
				}
			}
			prev_position = pos;
			}

			if(CigarCluster.recordDepthByPosition){
				map.addToEntry(pos, 1);
				maps[src_index].addToEntry(pos, 1);
				if(!match){
					errors[src_index].addToEntry(pos, 1);
				}
			}
		}

		public String toString() {
			return this.breaks.toString();
		}


		
		
		Integer getDepth(Integer i) {
			if(map==null) return 0;
			return this.map.getDepth(i);//this.map.containsKey(i) ?  map.get(i) :  0;
		}
		
		
		
		void getDepthSt(Integer i, int[] row, int start_pos, boolean match) {
			//StringBuffer sb = new StringBuffer();
			if(maps==null) return ;
			for(int src_index=0; src_index<maps.length; src_index++){
				row[start_pos+src_index] = match?  this.maps[src_index].getDepth(i) : this.errors[src_index].getDepth(i) ;
			}
		}
		
		String getTotDepthSt(boolean match) {
			if(maps==null) return "NA";
			StringBuffer sb = new StringBuffer();
			for(int src_index=0; src_index<maps.length; src_index++){
				if(src_index>0)sb.append("\t");
				sb.append(match ? this.maps[src_index].valsum() : this.errors[src_index].valsum());
			}
			return sb.toString();
		}
		public int readCountSum() {
			return readCountSum;
		}
		
		String getError(int src_index){
			if(errors==null) return "NA";
			double err = (double)this.errors[src_index].valsum()/(double)this.maps[src_index].valsum();
			if(err<-1e-5 || err> 1.0001) throw new NumberFormatException(" error is outside range of 0 1"+errors[src_index].valsum()+" "+maps[src_index].valsum());
			String st = this.maps[src_index].valsum()==0 ?  "NaN" :  String.format("%5.3g", err);
			return st;
		}
		
		String getErrorRatioSt() {
			if(errors==null) return "NA";
			StringBuffer sb = new StringBuffer();
			for(int src_index=0; src_index<maps.length; src_index++){
				if(src_index>0)sb.append("\t");
				sb.append(this.maps[src_index].valsum()==0 ?  "NaN" :  String.format("%5.3g", (double)this.errors[src_index].valsum()/(double)this.maps[src_index].valsum()));
			}
			return sb.toString();
		}
		
		
		int[][] exons;
		
		void addZeros(int seqlen){
			if(map==null) return;
			List<Integer> keys = this.map.keys();
			if(start> 1){
				map.addToEntry(start-1, 0);
			}
			for(int i=1; i<keys.size(); i++){
				if(keys.get(i)-keys.get(i-1)> 10){
					map.addToEntry(keys.get(i)-1, 0);
					map.addToEntry(keys.get(i-1)+1, 0);
				}
			}
			if(end<seqlen){
				map.addToEntry(end+1, 0);
			}
			
			
		}
		
		public int[][] getClusterDepth(int num_sources, Sequence ref) {
			//numPos =0;
			totLen =0;
			int numcols = 1 + num_sources * 2 + 1; 
		//	boolean prev0 = start>1;
			//boolean printPrev = false;
			List<Integer> keys = this.map.keys();
			int[][] matr = new int[keys.size()][numcols];
			int keysize = keys.size();
			
			for(int i=0; i<keysize; i++){
				Integer pos = keys.get(i);
				int[] row = matr[i];
				row[0] = pos;
				row[1] = (int) ref.getBase(pos-1); // because sequence is in 0 in index
				//System.err.println(ref.charAt(pos-1)+" -> " + row[1]);
				//A =, C = 1, G = 2. T = 3
				getDepthSt(pos, row,2, true);
				getDepthSt(pos, row,2 + num_sources, false);
				
			}
			return matr;
			//TranscriptUtils.printedLines[index]+=keys.size();
		}
		public int[][] getExons( double threshPerc, int numsteps) {
			List<Integer> start1 = new ArrayList<Integer>();
			List<Integer> end1 = new ArrayList<Integer>();
		
			double thresh = (double) readCountSum*threshPerc;
			boolean in =false;
			if(exons!=null) return exons;
			outer: for(int i=start; i<=end; i++) {
				double dep = getDepth(i);
				if(!in && dep>=thresh) {
					for(int j = 1; j<numsteps && i+j < end; j++) {
						if(getDepth(i+j)<thresh) continue outer; // no longer jumping in
					}
					in = true; 
					start1.add(i);
				}
				if(in && dep<thresh) {
					for(int j = 1; j<numsteps && i+j < end; j++) {
						if(getDepth(i+j)>=thresh) continue outer; // no longer jumping out
					}
					in  = false;
					end1.add(i-1);
				}
			}
			if(end1.size() < start1.size()) end1.add(end);
			 exons = new int[end1.size()][];
			for(int i=0; i<end1.size(); i++) {
				exons[i] = new int[] {start1.get(i), end1.get(i)};
				totLen += end1.get(i) - start1.get(i)+1;
			}
			return exons;
		}
	
		
		
		int[] readCount; 
		private int readCountSum;
		//int numPos =-1;
		int totLen = -1;
		public Collection<Integer> span = new TreeSet<Integer>();
		
		
		/** if its going to be less than thresh we return zero */
		public double similarity(CigarCluster c1, double thresh) {
			//if(this.index !=index) return 0;
			int overlap = TranscriptUtils.overlap(c1.start,c1.end,start, end);
		    if(overlap<0) return 0;
		    else{
		    	double union = (double) TranscriptUtils.union(c1.start, c1.end, start, end, overlap) ;
		    	if(overlap / union < thresh) return 0;
		    }
		//	double sim = map100.similarity(c1.map100);//this.similarity(map100, c1.map100);
			//if(sim<thresh ) return 0;
			double sim = map.similarity( c1.map);
			//System.err.println(highRes+" "+sim);
			return sim;
		}
		
		
		
		public static double similarity(int[][] exons1, int[][] exons2 ){
			int overlap =0;
			int union =0;
			for(int i=0; i<exons1.length; i++){
				int st1 = exons1[i][0];
				int end1 = exons1[i][1];
				for(int j=0; j<exons2.length; j++){
					int st2  =  exons2[j][0];
					int end2 = exons2[j][1];
					int overl = TranscriptUtils.overlap(st1, end1,st2, end2);
					int unio = TranscriptUtils.union(st1, end1,st2, end2, overl);
					overlap+=overl;
					union+= unio;
				}
			}
			return (double) overlap/(double) union;
		}
		public double exonSimilarity(CigarCluster c1){
			return similarity(c1.exons, this.exons);
		}
		public double similarity(CigarCluster c1) {
			int overlap = TranscriptUtils.overlap(c1.start, c1.end,start, end);
		    if(overlap<=0) return 0;
		    int union = TranscriptUtils.union(c1.start, c1.end,start, end, overlap);
		    double sim = (double) overlap/(double) union;
	    	if(sim < 0.5) return 0;
	    	//double sim1 = map100.similarity(c1.map100);
	    	//if(sim1<0.5) return sim1;
	    	 return  map.similarity( c1.map);  
		}
		
		
	/*	static int sum(Map<Integer, Integer> source){
			return source.values().stream() .reduce(0, Integer::sum);
		}*/
		
		public CigarHash2 merge(CigarCluster c1, int num_sources, int src_index) {
			if(c1.start < start) start = c1.start;
			if(c1.end > end) end = c1.end;
			if(breakSt<0 || (c1.breakSt>=0 & c1.breakSt > breakSt)) breakSt = c1.breakSt;
			if(breakEnd<0 || (c1.breakEnd>=0 &  c1.breakEnd < breakEnd)) breakEnd = c1.breakEnd;
			if(breakSt2<0 || (c1.breakSt2>=0 & c1.breakSt2 > breakSt2)) breakSt2 = c1.breakSt2;
			if(breakEnd2<0 || (c1.breakEnd2>=0 &  c1.breakEnd2 < breakEnd2)) breakEnd2 = c1.breakEnd2;
			//int subID;
			this.span.addAll(c1.span);
			if(c1.all_breaks!=null){
				/*for(Iterator<CigarHash2> it = c1.all_breaks.keySet().iterator() ; it.hasNext();){
					CigarHash2 nxt = it.next();
					Integer count = this.all_breaks.get(nxt);
					all_breaks.put(nxt, count==null ? 1 : count+1);
				}*/
				throw new RuntimeException("!!");
			}
				
			CigarHash2 br = (CigarHash2) c1.breaks.clone(true);
			Count	count = this.all_breaks.get(br);
				if(count==null) {
					count = new Count(num_sources, src_index, all_breaks.size());
					all_breaks.put(br, count);
				}
				else count.increment(src_index);
			
			for(int i=0; i<this.readCount.length;i++) {
				readCount[i]+=c1.readCount[i];
			}
			this.readCountSum+=c1.readCountSum;
			//System.err.println(this.breaks.toString()+"\t"+c1.index+" "+readCountSum);
					
			if(this.recordDepthByPosition){
				map.merge( c1.map);
				for(int i=0; i<maps.length; i++){
					maps[i].merge( c1.maps[i]);
					errors[i].merge(c1.errors[i]);
				}
			}
			return br;
		}
	
		

		public int numBreaks() {
			return (int)Math.round(((double)this.breaks.size()-2.0)/2.0);
		}

		public String exonCount(){
			Set<Integer>s = new TreeSet<Integer>();
			
			for(Iterator<CigarHash2> it = this.all_breaks.keySet().iterator();it.hasNext();){
				CigarHash2 br = it.next();
				s.add((int) Math.round((double)br.size()/2.0));
			}
			StringBuffer sb = new StringBuffer();
			boolean first = true;
			for(Iterator<Integer> it = s.iterator(); it.hasNext(); ){
				if(first){
					first = false;
				}else{
					sb.append(",");
				}
				sb.append(it.next());
			}
			String st = sb.toString();
			return st;
		}
		
		
	}