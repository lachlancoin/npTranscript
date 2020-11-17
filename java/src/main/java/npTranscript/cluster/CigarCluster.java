package npTranscript.cluster;

import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
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
		public static int annotationToRemoveFromGFF = -1;
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
/** returns average break point position */
	public static class Count implements Comparable{
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
			if(true_breaks==null) return null;
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
		@Override
		public int compareTo(Object o) {
			return Integer.compare(this.sum(), ((Count)o).sum());
		}
		
		
	}
		
   	Map<CigarHash2, Count> all_breaks  ;
   	CigarHash2 breaks = new CigarHash2();
   	CigarHash breaks_hash = new CigarHash();
   	
   	//chr1    HAVANA  gene    11869   14409   .       +       .       
   	//ID=ENSG00000223972.5;gene_id=ENSG00000223972.5;gene_type=transcribed_unprocessed_pseudogene;gene_name=DDX11L1;level=2;havana_gene=OTTHUMG00000000961.2
private static void writeGFF1(List<Integer> breaks, PrintWriter pw,SequenceOutputStream os,  PrintWriter[] bedW, String chr, 
		String type,  String parent, String ID,  int start, int end, String type_nme, String secondKey, String geneID, char strand, Sequence seq,  int []counts){
	//String secondKey = this.breaks_hash.secondKey;
	 //char strand = '+';//secondKey.charAt(secondKey.length()-1);
	 int count_ref= CigarCluster.annotationToRemoveFromGFF>=0  ? counts[annotationToRemoveFromGFF]:0;
		if(count_ref>0) return;// do not print if count ref greater than zero;

	pw.print(chr);pw.print("\tnp\t"+type+"\t");
	 pw.print(start);pw.print("\t"); pw.print(end);
	 pw.print("\t.\t"); pw.print(strand);pw.print("\t.\t");
	 pw.print("ID="); pw.print(ID);
	 if(parent!=null) {
		 pw.print(";Parent=");pw.print(parent);
	 }
	 pw.print(";gene_id=");pw.print(geneID);
	 pw.print(";gene_type=");pw.print(type_nme);
	 pw.print(";type=ORF;");
	 pw.print("gene_name=");pw.print(secondKey.replace(';','_'  ));
	 int count =0 ;
	
	 for(int i=0 ; i<counts.length; i++){
		 count+=counts[i];
		
	 }
	 String tpmstr = count+"";//String.format("%5.3g", "tpm");
	 pw.print(";count=");pw.print(tpmstr);

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
		 pw.print(";gene_id=");pw.print(geneID);
		 pw.print(";type=ORF;");
		 pw.print("gene_name=");pw.print(secondKey.replace(';','_'  ));
		 pw.println();
	 }
	  Sequence seq1 = new Sequence(seq.alphabet(), sb.toString(), ID);
	  seq1.setDesc(chr+" "+breaks.toString()+" "+type_nme+" "+tpmstr);
	  try{
	  seq1.writeFasta(os);
	  }catch(Exception exc){
		  exc.printStackTrace();
	  }
	 int source_index =0;
		if(bedW!=null){
			for(int i=0; i<bedW.length; i++){
				if(counts[i]>0){
					printBed(bedW[i], chr, breaks, ID,  strand, source_index, geneID,   counts[i]);
				}
			}
		}
	  
	 }
}

public static  synchronized void printBed(PrintWriter bedW, String chrom, List<Integer> breaks, String transcript_id, char strand, int source, String gene_id,   int read_depth){//, int span, String span_str){
	if(bedW==null) return;
	int num_exons =(int) Math.floor( (double)  breaks.size()/2.0);
	int col_id = source % Outputs.col_len;

	int startPos = breaks.get(0)-1;
	int endPos = breaks.get(breaks.size()-1)-1;
	StringBuffer block_sizes = new StringBuffer();
	StringBuffer block_start = new StringBuffer();
	String comma = "";
	
	for(int i=0; i<breaks.size(); i+=2){
		block_start.append(comma+(breaks.get(i)-1-startPos));
		block_sizes.append(comma+(breaks.get(i+1)-breaks.get(i)));
		if(i==0) comma=",";
	}
	bedW.println(chrom+"\t"+startPos+"\t"+endPos+"\t"+transcript_id+"."+gene_id+"\t"+read_depth+"\t"+strand+"\t"+startPos+"\t"+endPos+"\t"
			+Outputs.col_str[col_id]+"\t"+num_exons+"\t"+block_sizes.toString()+"\t"+block_start.toString());

}

 public void writeGFF(PrintWriter pw, SequenceOutputStream os, PrintWriter[] bedW, String chr, double  iso_thresh, String type_nme, Sequence seq){
	// String secondKey =  this.breaks_hash.secondKey;
		if(this.readCountSum< Outputs.gffThreshGene) return;
//if(!type_nme.equals("5_3")) return;
	//double tpm = (double) this.readCountSum() / (double) mapped_read_count;
	writeGFF1(null, pw, os,null, chr, "gene",  null, this.breaks_hash.secondKey, this.start, this.end, type_nme, this.breaks_hash.secondKey, 
			this.id, strand, seq,this.readCount);
	if(all_breaks.size()==0) throw new RuntimeException("no transcripts");
	List<Count> counts= new ArrayList<Count>(all_breaks.values());
	Collections.sort(counts);
	int max_cnt = counts.get(counts.size()-1).sum();
	int min_cnt = counts.get(0).sum();
	//System.err.println(min_cnt+" "+max_cnt+" "+counts.size());
		 if(min_cnt >max_cnt) throw new RuntimeException();
	
		 //double thresh1 = Math.min(Outputs.gffThresh, b)
		inner: for(int i=counts.size()-1; i>=counts.size()-Outputs.maxTranscriptsPerGeneInGFF; i--) {
			Count br_next = counts.get(i);
			String gene_id =Outputs.maxTranscriptsPerGeneInGFF==1 ?  this.breaks_hash.secondKey:this.breaks_hash.secondKey +".t"+i;
			if(br_next.sum() < Outputs.gffThreshTranscript) break inner;
			
			List<Integer> br_ = br_next.getBreaks();
				writeGFF1(br_, pw, os ,bedW, chr, "transcript", 
						this.id, gene_id,  br_.get(0),
						br_.get(br_.size()-1), type_nme, this.breaks_hash.secondKey, this.id, strand,seq,  br_next.count);
				
				
			
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
				if(CigarCluster.recordStartEnd){
					this.mapStart = new SparseArrayVector(num_sources);
					this.mapEnd = new SparseArrayVector(num_sources);
				}
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
		public void setStartEnd(int startPos, int endPos, int src) {
			if(CigarCluster.recordStartEnd){
				this.mapStart.addToEntry(start, src,1);
				this.mapEnd.addToEntry(end, src, 1);
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
			this.breaks.addAllR(c1.breaks,0);
			all_breaks = new HashMap<CigarHash2, Count>();
			//if(all_breaks.size()>0) throw new RuntimeException("should be zero");
			Count cnt =new Count(num_sources, source_index,  0);
			if(Outputs.writeGFF || Outputs.writeIsoforms) cnt.addBreaks(c1.breaks);
			this.all_breaks.put(IdentityProfile1.includeStart  ? breaks : breaks.clone(true, 1),cnt);
			this.breaks_hash.setSecondKey(c1.breaks_hash.secondKey);
			
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
				if(recordStartEnd){
						mapStart.merge( c1.mapStart);
						mapEnd.merge(c1.mapEnd);
				}
			}
		}

		 final SparseVector map;// = new SparseVector(); //coverage at high res
		final  SparseVector[] maps, errors;
		SparseArrayVector mapStart, mapEnd;
		
		final private char strand;
	   public void clear(){
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
				if(this.recordStartEnd){
					this.mapStart.clear();
					this.mapEnd.clear();
				}
			}
			Arrays.fill(readCount, 0);
			
			readCountSum=1;
			
			this.prev_position = 0;
			this.break_point_cluster = -1;
			first = true;
			this.breaks.clear();			
			if(all_breaks!=null)this.all_breaks.clear();
			this.breaks_hash.clear();
	   }

		
	
public static boolean recordStartEnd = false;
	   
		public void add(int pos, int src_index, boolean match) {
			boolean break_p = pos-prev_position>IdentityProfile1.break_thresh;
			if(CigarCluster.recordDepthByPosition){
				
				map.addToEntry(pos, 1);
				maps[src_index].addToEntry(pos, 1);
				if(!match){
					errors[src_index].addToEntry(pos, 1);
				}
				if(recordStartEnd && break_p ){
					
					mapStart.addToEntry(pos,src_index, 1);
					mapEnd.addToEntry(prev_position,src_index, 1);

				}
			}
			
			if(match){
			//System.err.println(prev_position + " "+pos);
			if(first){
				first = false; 
				breaks.add(pos);
			} else{	
				if(break_p){
					this.breaks.add(prev_position);
					this.breaks.add(pos);
				}
			}
			prev_position = pos;
			}

			
		}

		public String toString() {
			return this.breaks.toString();
		}


		
		
		Integer getDepth(Integer i) {
			if(map==null) return 0;
			return this.map.getDepth(i);//this.map.containsKey(i) ?  map.get(i) :  0;
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
		
		
		String getTotDepthSt(boolean match) {
			if(maps==null) return "NA";
			StringBuffer sb = new StringBuffer();
			for(int src_index=0; src_index<maps.length; src_index++){
				if(src_index>0)sb.append("\t");
				sb.append(String.format("%5.3g", match ? this.maps[src_index].valsum() : this.errors[src_index].valsum()).trim());
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
		
		/*void addZeros(int seqlen){
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
			
			
		}*/
		
		public void getClusterDepth(int[][] matr, List<Integer> keys,   int[] col_inds, int offset) {
			totLen =0;
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
			
			CigarHash2 br = (CigarHash2) c1.breaks.clone(true,IdentityProfile1.includeStart ? 0 : 1);
			Count	count = this.all_breaks.get(br);
				if(count==null) {
					count = new Count(num_sources, src_index, all_breaks.size());
				
					all_breaks.put(br, count);
				}
				else{
				count.increment(src_index);
				}
			if(Outputs.writeGFF || Outputs.writeIsoforms) count.addBreaks(c1.breaks);
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
				if(recordStartEnd){
					
						mapStart.merge( c1.mapStart);
						mapEnd.merge(c1.mapEnd);
				}
			}
			return br;
		}
	
		

		public int numIsoforms() {
			return all_breaks==null ? 0 : (int)Math.round(((double)this.all_breaks.size()-2.0)/2.0);
		}
		public List<Integer> getBreaks() {
			int max=0;
			Count maxv = null;
			for(Iterator<Count> it = this.all_breaks.values().iterator();it.hasNext();){
				Count br = it.next();
				if(br.sum()>max){
					max = br.sum();
					maxv = br;
				}
			}
			if(maxv!=null){
				return maxv.getBreaks();
			}
			else return null;
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