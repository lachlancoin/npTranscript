package npTranscript.cluster;

import java.util.List;

public class Maps{
	final SparseVector map;// = new SparseVector(); //coverage at high res
	 final  SparseVector[] maps, errors;
	final SparseArrayVector mapStart, mapEnd;
	public void clear(){
		map.clear();
		this.totalDepth=1;
		for(int i=0; i<maps.length; i++){
			maps[i].clear();
			errors[i].clear();
		}
			this.mapStart.clear();
			this.mapEnd.clear();
	}
	int totalDepth=1;
	public int totalDepth() {
		// TODO Auto-generated method stub
		return totalDepth;
	}
	public void getClusterDepth(int[][] matr, List<Integer> keys,   int[] col_inds, int offset) {
		//totLen =0;
		 int keysize = keys.size();
		
		for(int i=0; i<keysize; i++){
			Integer pos = keys.get(i);
			int[] row = matr[i];
			row[0] = pos;
			getDepthSt(pos, row, col_inds, offset);
		}
	}
	public static void getClusterDepthStartEnd(int[][] matr, List<Integer> keys,   int offset, SparseArrayVector mapS) {
		 int keysize = keys.size();
		int len = mapS.len;
		for(int i=0; i<keysize; i++){
			Integer pos = keys.get(i);
			int[] row = matr[i];
			row[0] = pos;
			System.arraycopy(mapS.get(pos), 0, row, 1, len);
		}
	}
	
	void getDepthSt(Integer i, int[] row,   int[] col_inds, int offset) {
		//StringBuffer sb = new StringBuffer();
		if(maps==null) return ;
		for(int src_index=0; src_index<maps.length; src_index++){
			int i1 =offset+2*(col_inds[src_index]-offset );
			row[i1] =  this.maps[src_index].getDepth(i) ;
			row[i1+1]	=	this.errors[src_index].getDepth(i) ;
		}
	}
	
	
	public void merge(Maps c1) {
		//System.err.println("merge "+this.)
		this.totalDepth= this.totalDepth+c1.totalDepth();
	//	System.err.println("total depth "+totalDepth);
		map.merge( c1.map);
		for(int i=0; i<maps.length; i++){
			maps[i].merge( c1.maps[i]);
			errors[i].merge(c1.errors[i]);
			if(errors[i].valsum()<0) throw new NumberFormatException("shoud not be negative");
			if(errors[i].valsum()>maps[i].valsum()) throw new NumberFormatException("shoud not be greater");
		}
				mapStart.merge( c1.mapStart);
				mapEnd.merge(c1.mapEnd);
		
	}
	public void setStartEnd(int startPos, int endPos, int src) {
			this.mapStart.addToEntry(startPos, src,1);
			this.mapEnd.addToEntry(endPos, src, 1);
		
	}
	
	Maps(Maps maps){
		this(maps.num_sources);
		this.map.transferFrom(maps.map);
		this.mapStart.transferFrom(maps.mapStart);
		this.mapEnd.transferFrom(maps.mapEnd);
		for(int i=0; i<num_sources; i++){
			this.maps[i].transferFrom(maps.maps[i]);
			this.errors[i].transferFrom(maps.errors[i]);
		}
		
	}
	final int num_sources;
	Maps(int num_sources){
		this.num_sources = num_sources;
		map = new SparseVector();
		this.maps = new SparseVector[num_sources];
		this.errors= new SparseVector[num_sources];
		this.mapStart = new SparseArrayVector(num_sources);
		this.mapEnd = new SparseArrayVector(num_sources);
		for(int i=0; i<maps.length; i++){
			maps[i] = new SparseVector();
			errors[i] = new SparseVector();
		}
	}
	Integer getDepth(Integer i) {
		if(map==null) return 0;
		return this.map.getDepth(i);//this.map.containsKey(i) ?  map.get(i) :  0;
	}
	String getTotDepthSt(boolean match) {
		if(maps==null) return "NA";
		StringBuffer sb = new StringBuffer();
		for(int src_index=0; src_index<maps.length; src_index++){
			if(src_index>0)sb.append("\t");
			sb.append(String.format("%5.3g", match ? this.maps[src_index].valsum() : this.errors[src_index].valsum()).trim());
		}
		return sb.toString();
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
	public void addToEntry(int src_index, int pos, int i, boolean match, boolean break_p, int prev_position) {
		map.addToEntry(pos, 1);
		maps[src_index].addToEntry(pos, 1);
		if(!match){
			errors[src_index].addToEntry(pos, 1);
		}
		if(break_p && prev_position >0 ){
			
			mapStart.addToEntry(pos,src_index, 1);
			mapEnd.addToEntry(prev_position,src_index, 1);

		}
		
	}

	
	
	
}