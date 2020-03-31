package npTranscript.cluster;

import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.List;

import ch.systemsx.cisd.hdf5.IHDF5SimpleWriter;

/**
 * @author Lachlan Coin
 *
 */

public class CigarCluster  {
		//final int index;
		
		
		static int round2 = 100;
		int breakSt = -1;
		int breakEnd = -1;
		
		final int id;

		 int start=0;
		 int end=0;
		int prev_position =0;
		int break_point_cluster = -1;		
		boolean first = true;
		
	
		public void addReadCount(int source_index) {
			readCount[source_index]++;
			this.readCountSum++;
			
		}
		
		
   	CigarHash breaks = new CigarHash() ;
 
   	
   	

	
		
		public void setBreaks(List<Integer> breaks){
			this.breaks.addAll(breaks);

		}
		public CigarCluster(int id,  int num_sources){
			this.id = id;
			this.readCount = new int[num_sources];
			this.maps = new SparseVector[num_sources];
			this.errors= new SparseVector[num_sources];
			for(int i=0; i<maps.length; i++){
				maps[i] = new SparseVector();
				errors[i] = new SparseVector();
				
			}
		}
		
		public CigarCluster(int id,  int num_sources, CigarCluster c1, int source_index) {
			this(id, num_sources);
			this.breakSt = c1.breakSt;
			this.breakEnd = c1.breakEnd;
			setBreaks(c1.breaks);
			addReadCount(source_index);
			start = c1.start;
			end = c1.end;
			
			map.merge( c1.map);
			for(int i=0; i<maps.length; i++){
				maps[i].merge( c1.maps[i]);
				errors[i].merge(c1.errors[i]);
			}
		}

		private SparseVector map = new SparseVector(); //coverate at high res
	
		final private SparseVector[] maps, errors;
		public void clear(int source_index) {
			map.clear();
			this.breakSt=-1;
			this.breakEnd = -1;
			for(int i=0; i<maps.length; i++){
				maps[i].clear();
				errors[i].clear();
			}
			Arrays.fill(readCount, 0);
			readCount[source_index]=1;
			readCountSum=1;
			
			this.prev_position = 0;
			this.break_point_cluster = -1;
			first = true;
			this.breaks.clear();
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
			
			map.addToEntry(pos, 1);
		
			//int round1 = (int) Math.floor((double)pos/round2);
		//	map100.addToEntry(round(pos,round2),1);
		
			maps[src_index].addToEntry(pos, 1);
			
			if(!match){
				errors[src_index].addToEntry(pos, 1);
			}
		}

		public String toString() {
			return this.breaks.toString();
		}
/*
		public Iterator<Integer> keys() {
			return map100.keyIt();
		}
*/
		/*public Iterator<Integer> tailKeys(int st) {
			return map100.tailKeys(st);
		}*/

		
		
		Integer getDepth(Integer i) {
			return this.map.getDepth(i);//this.map.containsKey(i) ?  map.get(i) :  0;
		}
		
		void getDepthSt(Integer i, int[] row, int start_pos, boolean match) {
			//StringBuffer sb = new StringBuffer();
			for(int src_index=0; src_index<maps.length; src_index++){
				row[start_pos+src_index] = match?  this.maps[src_index].getDepth(i) : this.errors[src_index].getDepth(i) ;
			}
		}
		
		String getTotDepthSt(boolean match) {
			StringBuffer sb = new StringBuffer();
			for(int src_index=0; src_index<maps.length; src_index++){
				if(src_index>0)sb.append("\t");
				sb.append(match ? this.maps[src_index].valsum : this.errors[src_index].valsum);
			}
			return sb.toString();
		}
		
		String getError(int src_index){
			String st = this.maps[src_index].valsum==0 ?  "NaN" :  String.format("%5.3g", (double)this.errors[src_index].valsum/(double)this.maps[src_index].valsum);
			return st;
		}
		
		String getErrorRatioSt() {
			StringBuffer sb = new StringBuffer();
			for(int src_index=0; src_index<maps.length; src_index++){
				if(src_index>0)sb.append("\t");
				sb.append(this.maps[src_index].valsum==0 ?  "NaN" :  String.format("%5.3g", (double)this.errors[src_index].valsum/(double)this.maps[src_index].valsum));
			}
			return sb.toString();
		}
		
		
		int[][] exons;
		
		void addZeros(){
			boolean prev0 = start>1;
			boolean printPrev = false;
			//numPos =0;
		
			for(int i=this.start; i<this.end; i++) {
				int depth_i = getDepth(i);
				if(depth_i>0){
					if(prev0 && !printPrev){
					//	numPos++;
						this.map.addZero(i-1);
					}
					//numPos++;
					prev0 = false;
					printPrev = true;
				}else if(!prev0){
				//	numPos++;
					printPrev = true;
					prev0=true;
				}else{
					printPrev = false;
					prev0 = true;
				}
			}
			//return numPos;
		}
		
		public int[][] getClusterDepth(int num_sources) {
			//numPos =0;
			totLen =0;
			int numcols = 1 + num_sources * 2;
		//	boolean prev0 = start>1;
			//boolean printPrev = false;
			List<Integer> keys = this.map.keys();
			int[][] matr = new int[keys.size()][numcols];
			int keysize = keys.size();
			
			for(int i=0; i<keysize; i++){
				Integer pos = keys.get(i);
				int[] row = matr[i];
				row[0] = pos;
				getDepthSt(i, row,1, true);
				getDepthSt(i, row,1 + num_sources, false);
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
	
		
		
		int[] readCount; int readCountSum;
		//int numPos =-1;
		int totLen = -1;
		
		
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
		
		public void merge(CigarCluster c1) {
			if(c1.start < start) start = c1.start;
			if(c1.end > end) end = c1.end;
			if(breakSt<0 || (c1.breakSt>=0 & c1.breakSt < breakSt)) breakSt = c1.breakSt;
			if(breakEnd<0 || (c1.breakEnd>=0 &  c1.breakEnd > breakEnd)) breakEnd = c1.breakEnd;
			for(int i=0; i<this.readCount.length;i++) {
				readCount[i]+=c1.readCount[i];
			}
			this.readCountSum+=c1.readCountSum;
			//System.err.println(this.breaks.toString()+"\t"+c1.index+" "+readCountSum);
					
			 map.merge( c1.map);
		//	int sum2 = map100.merge(c1.map100);
			for(int i=0; i<maps.length; i++){
				maps[i].merge( c1.maps[i]);
				errors[i].merge(c1.errors[i]);
			}
			//if(sum1!=sum2){
				//throw new RuntimeException("maps not concordant");
			//}
		}
		
	}