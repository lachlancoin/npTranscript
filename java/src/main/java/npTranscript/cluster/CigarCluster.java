package npTranscript.cluster;

import java.io.PrintWriter;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;

import japsa.seq.Sequence;
import japsa.seq.SequenceOutputStream;


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

		
		}
		private int id;
		private int[] count;
		
		private double[]  true_breaks; //in 100kb
		private static double divisor = 1e5;
		public void addBreaks(List<Integer>breaks){
			if(true_breaks==null) true_breaks =new double[breaks.size()];
			if(true_breaks.length!=breaks.size()){
				System.err.println("warning breaks have different lengths, so not updating");
			}else{
				for(int i=0; i<true_breaks.length; i++){
					true_breaks[i] += breaks.get(i)/divisor;
				}
			}
		}
		public List<Integer> getBreaks(){
			Integer[] res = new Integer[true_breaks.length];
			int sum = this.sum();
			double mult = (divisor/sum);
			for(int i=0; i<res.length; i++){
				res[i] = (int) Math.round(true_breaks[i]*mult);
			}
			return Arrays.asList(res);
		}
		
		public void increment(int source_index) {
			this.count[source_index] = this.count[source_index]+1;
			
		}
		public int[] count(){
			return count;
		}
		public int id(){
			return id;
		}
		public int sum() {
			int sum = 0;
			for(int i=0; i<count.length; i++){
				sum+=count[i];
			}
			return sum;
		}
		
	}
		
   	Map<CigarHash2, Count> all_breaks  ;
   	CigarHash2 breaks = new CigarHash2();
   	CigarHash breaks_hash = new CigarHash();
   	
   	//chr1    HAVANA  gene    11869   14409   .       +       .       
   	//ID=ENSG00000223972.5;gene_id=ENSG00000223972.5;gene_type=transcribed_unprocessed_pseudogene;gene_name=DDX11L1;level=2;havana_gene=OTTHUMG00000000961.2
private static void writeGFF1(List<Integer> breaks, PrintWriter pw,SequenceOutputStream os,  String chr, 
		String type,  String parent, String ID,  int start, int end, String type_nme, String secondKey, String geneID, char strand, Sequence seq){
	//String secondKey = this.breaks_hash.secondKey;
	 //char strand = '+';//secondKey.charAt(secondKey.length()-1);
	 pw.print(chr);pw.print("\tnp\t"+type+"\t");
	 pw.print(start);pw.print("\t"); pw.print(end);
	 pw.print("\t.\t"); pw.print(strand);pw.print("\t.\t");
	 pw.print("ID="); pw.print(ID);
	 pw.print(";gene_id=");pw.print(geneID);
	 pw.print(";gene_type=");pw.print(type_nme);
	 pw.print(";type=ORF;");
	 pw.print("gene_name=");pw.print(secondKey);
	 if(parent!=null) {
		 pw.print(";Parent=");pw.print(parent);
	 }
	 pw.println();
	 if(breaks!=null){
		StringBuffer sb = new StringBuffer(); 
	for(int i=0; i<breaks.size();i+=2){
		int starti = breaks.get(i);
		int endi = breaks.get(i+1);
		sb.append( seq.subSequence(starti, endi).toString());
		 pw.print(chr);pw.print("\tnp\texon\t");
		 pw.print(starti);pw.print("\t"); pw.print(endi);
		 pw.print("\t.\t"); pw.print(strand);pw.print("\t.\t");
		 pw.print("ID="); pw.print(ID);pw.print(".e."+i);
		 pw.print(";Parent=");pw.print(ID);
		 pw.println();
	 }
	  Sequence seq1 = new Sequence(seq.alphabet(), sb.toString(), ID);
	  seq1.setDesc(chr+" "+breaks.toString()+" "+secondKey+" "+type_nme);
	  try{
	  seq1.writeFasta(os);
	  }catch(Exception exc){
		  exc.printStackTrace();
	  }
	 }
}


 public void writeGFF(PrintWriter pw, SequenceOutputStream os, String chr, double  iso_thresh, String type_nme, Sequence seq){
	// String secondKey =  this.breaks_hash.secondKey;
		if(this.readCountSum< Outputs.gffThresh) return;
//if(!type_nme.equals("5_3")) return;
	this.writeGFF1(null, pw, os, chr, "gene",  null, this.id, this.start, this.end, type_nme, this.breaks_hash.secondKey, this.id, strand, seq);
	Iterator<Count> it = this.all_breaks.values().iterator();
	double max = 0;

	while(it.hasNext()){
		int sum = it.next().sum();
		if(sum>max){
			max = sum;
		}
	}
	double minv = iso_thresh * max-0.001;
//	pw.print(breaks.get(0));pw.print("\t"); pw.print(breaks.get(breaks.size()));
		 Iterator<Count> it1 = this.all_breaks.values().iterator();
		for(int i=0; it1.hasNext();) {
			Count br_next = it1.next();
			List<Integer> br_ = br_next.getBreaks();
			if(br_next.sum()>=minv){
				writeGFF1(br_, pw, os ,chr, "transcript", this.id, this.id+".t"+i,  br_.get(0), br_.get(br_.size()-1), type_nme, this.breaks_hash.secondKey, this.id, strand,seq);
			}
		}
 }
	
		
		/*public void setBreaks(CigarHash2 breaks){
			this.breaks.addAll(breaks);
		}*/

		public CigarCluster(String id,  int num_sources, char strand){
			this.id = id;
			this.strand = strand;
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
		
		public CigarCluster(String id,  int num_sources, CigarCluster c1, int source_index, char strand) throws NumberFormatException{
			this(id, num_sources,strand);
			//this.strand = strand;
			this.span.addAll(c1.span);
			this.forward = c1.forward;
			this.breakSt = c1.breakSt;
			this.breakEnd = c1.breakEnd;
			this.breakSt2 = c1.breakSt2;
			this.breakEnd2 = c1.breakEnd2;
			this.breaks.addAllR(c1.breaks);
			all_breaks = new HashMap<CigarHash2, Count>();
			//if(all_breaks.size()>0) throw new RuntimeException("should be zero");
			Count cnt =new Count(num_sources, source_index,  0);
			if(Outputs.writeGFF) cnt.addBreaks(c1.breaks);
			this.all_breaks.put(breaks,cnt);
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
final private char strand;
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
		
		
		//int[][] exons;
		
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
			int seqlen = ref.length();
			int numcols = 1 + num_sources * 2 + 1; 
		//	boolean prev0 = start>1;
			//boolean printPrev = false;
			List<Integer> keys = this.map.keys(IdentityProfile1.writeCoverageDepthThresh);
			int[][] matr = new int[keys.size()][numcols];
			 int keysize = keys.size();
			
			for(int i=0; i<keysize; i++){
				Integer pos = keys.get(i);
				int[] row = matr[i];
				row[0] = pos;
				row[1] = pos <= ref.length() ? (int) ref.getBase(pos-1) : -1; // because sequence is in 0 in index
				//System.err.println(ref.charAt(pos-1)+" -> " + row[1]);
				//A =, C = 1, G = 2. T = 3
				getDepthSt(pos, row,2, true);
				getDepthSt(pos, row,2 + num_sources, false);
			}
			return matr;
			//TranscriptUtils.printedLines[index]+=keys.size();
		}
		
		
		
		int[] readCount; 
		private int readCountSum;
		//int numPos =-1;
		int totLen = -1;
		public Collection<Integer> span = new TreeSet<Integer>();
		
		
		
		
		
		
		
		
		
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
			//this.strand = c1.strand;
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
				else{
				count.increment(src_index);
				}
			if(Outputs.writeGFF) count.addBreaks(c1.breaks);
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