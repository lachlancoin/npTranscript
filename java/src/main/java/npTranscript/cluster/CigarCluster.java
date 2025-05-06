package npTranscript.cluster;

import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;
import java.util.TreeSet;
import java.util.stream.IntStream;

import japsa.seq.Sequence;
import japsa.seq.SequenceOutputStream;
import npTranscript.cluster.Outputs.HDFObj;


/**
 * @author Lachlan Coin
 *
 */

public class CigarCluster  {
		//final int index;
	
		
		public static boolean singleGFF = true;
	
		static int round2 = 100;
		public static boolean recordDepthByPosition = false;
		//public static int annotationToRemoveFromGFF = -1;
	
		/*int breakSt = -1;
		int breakEnd = -1;
		int breakSt2 = -1;
		int breakEnd2 = -1; */
		
		 Boolean forward = null;
		
		private final String id;
		
		final public String chrom;
		
		public String id(){
			return id;
		}
		/*public String id(int chrom_index){
			return chrom_index+"."+id;
		}*/

		 int startPos=0;
		 int endPos=0;
		int prev_position =-1; // if negative then indicates that it must be first
		int prev_ind = -1;
		int break_point_cluster = -1;		
	//	boolean first = true;
		boolean fusion = false;
		
	
		public void addReadCount(int source_index) {
			readCount[source_index]++;
			this.readCountSum++;
			
		}
/** returns average break point position */
	public static class Count implements Comparable{
		public Count(int[] count, String id){
			this.count = count;
			this.id = id;
		}
		public Count(int num_sources, int src_index, String id) {
			this.count = new int[num_sources];
			count[src_index]=1;
			this.id = id;
		
		}
		public Count(String id, HDFObj obj) {
			this.true_breaks = obj.br;
			this.count = obj.cnts;
			this.break_sum =  IntStream.of(this.count).sum(); 
			this.id = id;
		}
		private String id;
		private int[] count;
		
		public float[]  true_breaks; //in 1kb
		public static double divisor = 1e3f;
		
		int break_sum=0;
		
		
		//public static boolean baseBreakPointsOnFirst = true;
		
		public void addBreaks(List<Integer>breaks, int src_index){
		//	if(!Outputs1.firstIsTranscriptome || src_index==0 || this.count[0]==0){
				break_sum+=1;
				if(true_breaks==null) true_breaks =new float[breaks.size()];
				if(true_breaks.length!=breaks.size()){
					System.err.println("warning breaks have different lengths, so not updating");
				}else{
					for(int i=0; i<true_breaks.length; i++){
						true_breaks[i] += breaks.get(i)/divisor;
					}
				}
		//	}
		}
		
		
		public List<Integer> getBreaks(){
			if(true_breaks==null) return null;
			Integer[] res = new Integer[true_breaks.length];
			int sum = this.break_sum;//this.sum();
			double mult = (divisor/(double) sum);
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
		public String id(){
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
		
   	public Map<CigarHash2, Count> all_breaks  ;
   	CigarHash2 breaks = new CigarHash2(false);
   	List<Integer> start_positions = new ArrayList<Integer>();
   	CigarHash breaks_hash = new CigarHash();
   	
   	//chr1    HAVANA  gene    11869   14409   .       +       .       
   	//ID=ENSG00000223972.5;gene_id=ENSG00000223972.5;gene_type=transcribed_unprocessed_pseudogene;gene_name=DDX11L1;level=2;havana_gene=OTTHUMG00000000961.2
private static void writeGFF1(List<Integer> breaks,List<String> exons,  String key, PrintWriter pw,SequenceOutputStream os,  PrintWriter[] bedW, String chr, 
		String type,  String parent, String ID,  int start, int end, String type_nme, String secondKey, 
		String geneID, String strand, Sequence seq,  int []counts, int[] read_max){
	byte[] seqb =  seq==null ? null : seq.toBytes();
	boolean writeSeq = seqb!=null && seqb.length>0;
	 pw.print(chr);pw.print("\tnp\t"+type+"\t");
	 pw.print(start);pw.print("\t"); pw.print(end);
	 pw.print("\t.\t"); pw.print(strand);pw.print("\t.\t");
	 pw.print("ID="); pw.print(ID);
	 if(parent!=null) {
		 pw.print(";Parent=");pw.print(parent);
	 }
	 pw.print(";gene_id=");pw.print(parent==null ? ID: parent);
	 pw.print(";gene_type=");pw.print(type_nme);
	 pw.print(";type=ORF;");
	 if(key!=null){
		 pw.print(";key="+key+";"); 
	 }
	 pw.print("gene_name=");pw.print(secondKey.replace(';','_'  ));
	 int count =0 ;
	 pw.print(";count=");
	 for(int k=0; k<counts.length; k++){
		 count+=counts[k];
		 pw.print(counts[k]);
		 if(k<counts.length-1) pw.print(",");
	 }
	
	// String tpmstr = count+"";//String.format("%5.3g", "tpm");
	

	 pw.println();
	 if(breaks!=null){
		StringBuffer sb = new StringBuffer(); 
	for(int i=0; i<breaks.size();i+=2){
		int starti = breaks.get(i);
		int endi = breaks.get(i+1);
		if(writeSeq) sb.append( seq.subSequence(starti, endi).toString());
		 pw.print(chr);pw.print("\tnp\texon\t");
		 pw.print(starti);pw.print("\t"); pw.print(endi);
		 pw.print("\t.\t"); pw.print(strand);pw.print("\t.\t");
		 pw.print("ID="); pw.print(ID);pw.print(".e."+i);
		 pw.print(";Parent=");pw.print(ID);
		 pw.print(";gene_id=");pw.print(geneID);
		 pw.print(";type=ORF;");
		 pw.print("gene_name=");pw.print(secondKey.replace(';','_'  ));
		 pw.print(";count=");
		for(int k=0; k<counts.length; k++){
			 pw.print(counts[k]);
			 if(k<counts.length-1) pw.print(",");
		 }
		 pw.println();
	 }
	if(writeSeq){
	  Sequence seq1 = new Sequence(seq.alphabet(), sb.toString(), secondKey+"/"+ID);				
	  String read_count1 =  TranscriptUtils.getString(counts);

	  seq1.setDesc(chr+" "+CigarHash2.getString(breaks)+" "+CigarHash2.getString(exons)+" "+read_count1+" "+ IntStream.of(counts).sum() );//+" "+tpmstr);
	  try{
	  seq1.writeFasta(os);
	  }catch(Exception exc){
		  exc.printStackTrace();
	  }
	}
	 int source_index =0;
		if(bedW!=null){
			for(int i=0; i<bedW.length; i++){
				if(i<counts.length && counts[i]>0){
					printBed(bedW[i], chr, breaks, ID,  strand, i, geneID,   counts[i], read_max[i]);
				}
			}
		}
	  
	 }
}

public static  synchronized void printBed(PrintWriter bedW, String chrom, List<Integer> breaks, String transcript_id, 
		String strand, int source, String gene_id,   int read_depth, int read_max){//, int span, String span_str){
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
	int score =(int) Math.round(1000.0*((double)read_depth/(double)read_max));
	bedW.println(chrom+"\t"+startPos+"\t"+endPos+"\t"+transcript_id+"."+gene_id+"\t"+score+"\t"+strand+"\t"+startPos+"\t"+endPos+"\t"
			+Outputs.col_str[col_id]+"\t"+num_exons+"\t"+block_sizes.toString()+"\t"+block_start.toString());

}
static Comparator entryComparator = new Comparator<Entry<CigarHash2, Count>>(){
	@Override
	public int compare(Entry<CigarHash2, Count> o1, Entry<CigarHash2, Count> o2) {
		return o1.getValue().compareTo(o2.getValue());
	}
	
};

 public void writeGFF(PrintWriter[] pw, SequenceOutputStream os, PrintWriter[] bedW, String chr, double  iso_thresh, 
		 String type_nme, Sequence seq, Annotation annot){
	if(all_breaks.size()==0) throw new RuntimeException("no transcripts");
	List<Entry<CigarHash2,Count>> counts= new ArrayList<Entry<CigarHash2,Count>>(all_breaks.entrySet());
	Collections.sort(counts, entryComparator);
	
	int max_cnt = counts.get(counts.size()-1).getValue().sum();
	int min_cnt = counts.get(0).getValue().sum();
	boolean[] writeGene = new boolean[pw.length]; 
	boolean[] writeAny = new boolean[pw.length];
	int[] starts = new int[pw.length];
	int[] ends = new int[pw.length];
	int[] read_max = new int[bedW.length];
	for(int i=0; i<counts.size(); i++){
		int[] cnt=counts.get(i).getValue().count;
		for(int k=0; k<read_max.length; k++){
			read_max[k] = Math.max(read_max[k] ,cnt[k]);
 		}
	}
	Arrays.fill(starts, Integer.MAX_VALUE);;
	//System.err.println(min_cnt+" "+max_cnt+" "+counts.size());
	if(max_cnt >=1){
		 if(min_cnt >max_cnt) throw new RuntimeException();
		 
		 //double thresh1 = Math.min(Outputs.gffThresh, b)
		inner: for(int i=counts.size()-1; i>=counts.size()-Outputs1.maxTranscriptsPerGeneInGFF; i--) {
			if(i>=0){
			//	String keyv1 = counts.get(i).getKey().toString();
			//String keyv=	null;//counts.get(i).getKey().rescale().toString();
			
			Count br_next = counts.get(i).getValue();
			String keyv=br_next.id;
			int[] cnt = br_next.count;
			String gene_id =br_next.id();
			//int firstNonZero =num_sources-1;
			//boolean incl = false;
			Arrays.fill(writeGene, false);
			int thresh0  = Outputs1.gffThresh[0];
			for(int j=0; j<cnt.length; j++){
				if(cnt[j] >= thresh0) {
					//if any greater than thresh, then write all
					Arrays.fill(writeGene, true);
				}
			}
			if(Outputs1.firstIsTranscriptome){
				//then first depends only on first cnt
				int thresh1 = Outputs1.gffThresh[1];

				boolean someNonZero = false;
				for(int j=1; j<cnt.length; j++){
					if(cnt[j] >= thresh1) {
						someNonZero =true;
					}
				}
				writeGene[0] = cnt[0] >=1 && someNonZero;
			}
			
		

			List<Integer> br_ = br_next.getBreaks();
			List<String> exons = annot.matchExons(br_, this.chrom, forward);
			 for(int k=0; k<pw.length; k++){
				 if(writeGene[k]){
					 writeAny[k] = true;
					// System.err.println(k+this.id);
					int  start =  br_.get(0);;
					int end =  br_.get(br_.size()-1);
					if(start < starts[k]) starts[k] = start;
					if(end>ends[k]) ends[k] = end;
				writeGFF1(br_, exons, keyv , pw[k], os ,bedW, chr, "transcript", 
						this.id, gene_id, start,end
						, type_nme, this.breaks_hash.secondKey,gene_id,  strand,seq,  br_next.count, read_max);
				 }
			 }
			}
			
		}
		 for(int i=0; i<pw.length; i++){
			 if(writeAny[i]){
				 writeGFF1(null, null, null, pw[i], os,null, chr, "gene",  null, this.id, starts[i], ends[i], type_nme, this.breaks_hash.secondKey, 
							this.id, strand, seq,this.readCount, read_max);
			 }
		 }
	}
 }
 
 
		
		/*public void setBreaks(CigarHash2 breaks){
			this.breaks.addAll(breaks);
		}*/

		public CigarCluster(String chrom,String id,  int num_sources, String strand){
			this.id = id;
			this.chrom = chrom;
			this.strand = strand;
			this.readCount = new int[num_sources];
			if(recordDepthByPosition){
				this.map = new Maps(num_sources);
				
			}else{
				map = null;
			}
		}
		public void setStartEnd(int startPos, int endPos, int src) {
			if(map!=null){
				map.setStartEnd(startPos, endPos, src);
			}
			
		}
		
		
		public CigarCluster( String chrom, String strand, String key, String[] ids, Outputs.HDFObj[] obj){
			// /transcripts/U13369.1,NC_045512.2/296/37_40;291_296	
			//	 String[] add = addr[0].split("/");
			
				 this.chrom = chrom;// add[2];
				 this.id =chrom+"/"+strand+"/"+key;// add[2]+add[3];
				 this.strand = strand;
				// this.strand=key.charAt(key.length()-1);
				 int num_sources = obj[0].cnts.length;
				 this.readCount = new int[num_sources];
				 all_breaks = new HashMap<CigarHash2, Count>();	 
				 if(strand.startsWith("+")) forward = true;
				 else if(strand.startsWith("-")) forward= false;
				 this.breaks_hash.setSecondKey(id);
				 map = null;
				this.startPos = Integer.MAX_VALUE;
				this.endPos = -1;
				 for(int i=0; i<ids.length; i++){
					 Count cnt = new Count(ids[i],obj[i]);
					 int len_c = cnt.count.length;
					 if(len_c>num_sources){
						 num_sources=len_c;
						int[] rc1 = new int[len_c];
						System.arraycopy(this.readCount, 0, rc1, 0,this.readCount.length );
						this.readCount = rc1;
					 }
					 //System.err.println(len_c+" "+num_sources);
					for(int j=0; j<len_c; j++){
						this.readCount[j]+=cnt.count[j];
					}
					 String[] str = ids[i].split("[_,;]");
					 List<Integer> br = cnt.getBreaks();
					
					 this.startPos = Math.min(startPos, br.get(0));
					 this.endPos = Math.max(endPos, br.get(br.size()-1));
					 CigarHash2 ch = new CigarHash2(str);
					 this.all_breaks.put(ch,cnt);
					
				 }
				 
				 for(int j=0; j<num_sources; j++){
						this.readCountSum+=this.readCount[j];
					}
			}
				
		
		public CigarCluster(String chr, String id,  int num_sources, CigarCluster c1, int source_index, String strand,
				String subId) throws NumberFormatException{
			this(chr, id, num_sources,strand);
			//this.strand = strand;
			this.span.addAll(c1.span);
		//	if(c1.forward==null) throw new RuntimeException("!!");
			this.forward = c1.forward;
			/*this.breakSt = c1.breakSt;
			this.breakEnd = c1.breakEnd;
			this.breakSt2 = c1.breakSt2;
			this.breakEnd2 = c1.breakEnd2*/;
			all_breaks = new HashMap<CigarHash2, Count>();
			Count cnt =new Count(num_sources, source_index,  subId);
			this.breaks = c1.cloneBreaks();
			this.all_breaks.put(breaks,cnt);
			

			//if(all_breaks.size()>0) throw new RuntimeException("should be zero");
			if( Outputs.writeIsoforms) cnt.addBreaks(c1.breaks, source_index);
			//this.all_breaks.put(IdentityProfile1.includeStartEnd  || breaks.size()==2  ? breaks : breaks.clone(true, 1, breaks.size()-1),cnt);

			this.breaks_hash.setSecondKey(c1.breaks_hash.secondKey);
			
			addReadCount(source_index);
			startPos = c1.startPos;
			endPos = c1.endPos;
			if(map!=null){
				map.merge( c1.map);

				
			}
		}
final Maps map;
		 
		
	String strand;
	   public void clear(){
		   span.clear();
			this.forward = null;
			/*this.breakSt=-1;
			this.breakEnd = -1;
			this.breakSt2=-1;
			this.breakEnd2 = -1;*/
			if(map!=null){
				map.clear();
			}
			Arrays.fill(readCount, 0);
			
			readCountSum=1;
			
			this.prev_position = -1;
			this.break_point_cluster = -1;
		//	first = true;
			this.breaks.clear();	
			this.start_positions.clear();
			if(all_breaks!=null)this.all_breaks.clear();
			this.breaks_hash.clear();
	   }

	   public void addStart(int alignmentStart, int ref_ind) {
			// TODO Auto-generated method stub
			 this.start_positions.add(breaks.size());

			breaks.add(alignmentStart);
			this.prev_position=  alignmentStart;
			this.prev_ind = ref_ind;
		}
	
public static boolean recordStartEnd = false;
public int addBlock(int start, int len, boolean first, boolean last, int ref_ind, boolean negStrand) {
return addBlock(start, len, first, last, ref_ind,start-prev_position , negStrand);
}
	public int addBlock(int start, int len, boolean first, boolean last, int ref_ind, int gap, boolean negStrand) {
		int mult = 1;//negStrand ? -1 : 1;
		if(first) {
			this.start_positions.add(breaks.size());
		}
		int end = start + len-1;
		if(first && last){
			breaks.add(mult*start); breaks.add(mult*end);
		}
		else if(first){
			breaks.add(mult*start);
		}else if (last){
			breaks.add(mult*end);
		}else if(ref_ind !=prev_ind){
			breaks.add(mult*prev_position);
			breaks.add(mult*start);
		}else if(gap>IdentityProfile1.break_thresh){
			breaks.add(mult*prev_position);
			breaks.add(mult*start);
		}
		prev_position = end;
		this.prev_ind = ref_ind;
		return gap;
		
	}
		
		public void addEnd(int alignmentEnd, int ref_ind) {
			this.breaks.add(alignmentEnd);
			
		}
		public int add(int pos, int src_index, boolean match, int ref_ind) {
			return add(pos, src_index, match, ref_ind,pos-prev_position );
		}
		public int add(int pos, int src_index, boolean match, int ref_ind, int gap) {
		//	System.err.println(pos);
			boolean break_p =prev_position < 0 ||  gap>IdentityProfile1.break_thresh || ref_ind!=prev_ind;
			if(map!=null){
				map.addToEntry(src_index, pos, 1, match, break_p, prev_position);
			}
			
			if(match){
			//System.err.println(prev_position + " "+pos);
			
				if(break_p){
					if(prev_position>0) this.breaks.add(prev_position);
					this.breaks.add(pos);
				}
				prev_position = pos;
				this.prev_ind = ref_ind;
			}
			return gap;
			
		}

		public String toString() {
			return this.breaks.toString();
		}


		
		
		
		public int readCountSum() {
			return readCountSum;
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
		
		
		
		
		
		int[] readCount; 
		private int readCountSum;
		//int numPos =-1;
		int totLen = -1;
		public Collection<Integer> span = new TreeSet<Integer>();
		
		
		
		
		
		
		
		
		
	/*	static int sum(Map<Integer, Integer> source){
			return source.values().stream() .reduce(0, Integer::sum);
		}*/
		
		public CigarHash2 cloneBreaks(){
			if(IdentityProfile1.includeStartEnd ) return  (CigarHash2)  breaks.clone(true,0, breaks.size()) ;
			else if(forward!=null){ // we know strand
				if(forward){
					return breaks.clone(true,1,breaks.size());
				}else{
					return breaks.clone(true,0,breaks.size()-1);
				}
			}else{
				 return  (CigarHash2)  breaks.clone(true,0, breaks.size()) ;
			}
		}
		
		public CigarHash2 merge(CigarCluster c1, int num_sources, int src_index, String[] clusterID, String readId) {
			if(c1.startPos < startPos) startPos = c1.startPos;
			if(c1.endPos > endPos) endPos = c1.endPos;
			/*if(breakSt<0 || (c1.breakSt>=0 & c1.breakSt > breakSt)) breakSt = c1.breakSt;
			if(breakEnd<0 || (c1.breakEnd>=0 &  c1.breakEnd < breakEnd)) breakEnd = c1.breakEnd;
			if(breakSt2<0 || (c1.breakSt2>=0 & c1.breakSt2 > breakSt2)) breakSt2 = c1.breakSt2;
			if(breakEnd2<0 || (c1.breakEnd2>=0 &  c1.breakEnd2 < breakEnd2)) breakEnd2 = c1.breakEnd2;
			*///int subID;
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
			
			CigarHash2 br = c1.cloneBreaks();
				
			Count	count = this.all_breaks.get(br);
				if(count==null) {
					String id;
					if(Outputs1.firstIsTranscriptome && src_index==0){
						id = readId;
					}else{
						id =this.id+".t"+all_breaks.size();
					}
					count = new Count(num_sources, src_index, id);
				
					all_breaks.put(br, count);
				}
				else{
					if(Outputs1.firstIsTranscriptome && src_index==0){
						count.id = readId;
					}
				count.increment(src_index);
				}
				clusterID[1] = count.id()+"";
				//if(!clusterID[1].startsWith("R") && src_index==0){
				//	throw new RuntimeException("!!");
				//}
			if( Outputs.writeIsoforms){
				count.addBreaks(c1.breaks, src_index);
			}
			for(int i=0; i<this.readCount.length;i++) {
				readCount[i]+=c1.readCount[i];
			}
			this.readCountSum+=c1.readCountSum;
			//System.err.println(this.breaks.toString()+"\t"+c1.index+" "+readCountSum);
					
			if(map!=null){
				map.merge( c1.map);
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

		public String getBreakString(String strand) {
			// TODO Auto-generated method stub
			StringBuffer sb = new StringBuffer();
			boolean reverse = strand.charAt(0)=='+';
			if(reverse) {
				for(int i=this.start_positions.size()-1;i>=0; i--){
					int fromIndex = start_positions.get(i);
					
					//String st1;
					if(i < start_positions.size()-1){ 
						int toIndex = start_positions.get(i+1);
						sb.append(CigarHash2.getString(breaks.subList(fromIndex,  toIndex), strand.charAt(i))); 
						sb.append(";");	
					}else{
						int toIndex =  breaks.size();
						sb.append(CigarHash2.getString(breaks.subList(fromIndex,  toIndex), strand.charAt(i)));
					}
				}
			}else {
			for(int i=0; i<this.start_positions.size(); i++){
				int fromIndex = start_positions.get(i);
				
				//String st1;
				if(i < start_positions.size()-1){ 
					int toIndex = start_positions.get(i+1);
					sb.append(CigarHash2.getString(breaks.subList(fromIndex,  toIndex))); 
					sb.append(";");	
				}else{
					int toIndex =  breaks.size();
					sb.append(CigarHash2.getString(breaks.subList(fromIndex,  toIndex)));
				}
			}
			}
			return sb.toString();
		}
static List<String> empty_list = Arrays.asList(new String[0]);
		
		public void process1( Outputs1 o, Sequence genome, Annotation annot ){
		//	genomes==null ? null : genomes.get(chrom)
			CigarCluster cc = this;
			
			
			String read_count = TranscriptUtils.getString(cc.readCount);
			String chrom = cc.chrom;//seq.getName();
			
			Boolean forward = cc.forward;
		//	boolean hasLeaderBreak = false;//TranscriptUtils.coronavirus  ? (cc.breaks.size()>1 &&  annot.isLeader(cc.breaks.get(1)*CigarHash2.round)) : false;
			//geneNames.clear();
		//	int type_ind = 0;//annot.getTypeInd(cc.startPos, cc.endPos, forward);
			String type_nme ="NA";// annot.nmes[type_ind];
			//String geneNme = "NA";//annot.getString(cc.span, geneNames);
			if(Outputs1.writeGFF){
				cc.writeGFF(o.gffW, o.refOut, o.bedW,chrom,  Outputs.isoThresh, type_nme, genome, annot);
			}
		
		
			//		
				//		
			
			String depth_str ="";// "\t"+cc.getTotDepthSt(true)+"\t"+cc.getTotDepthSt(false);
			Iterator<Entry<CigarHash2, Count>> it = cc.all_breaks.entrySet().iterator();
			Set<String> genes = new TreeSet<String>();
			while(it.hasNext()){
				Entry<CigarHash2, Count> nxt = it.next();
				CigarHash2 h2 = nxt.getKey();
				Count cnt = nxt.getValue();
				List<Integer> breaks = cnt.getBreaks();
				 List<String> exons  = annot==null ? annot.empty_list : annot.matchExons(breaks, chrom, this.forward);
				genes.addAll(exons);
				//int exonCount =(int) Math.round((double) breaks.size()/2.0);
				int exonCount = exons.size();
			//QUESTION ABOUT HOW TO REPORT EXON STRING
				String exon_str =CigarHash2.getString(exons);//CigarHash2.getString1(exons," ;; ")
				String read_count1 =  TranscriptUtils.getString(cnt.count());
				String str = cc.id()+"/"+cnt.id()+"\t"+chrom+"\t"+breaks.get(0)+"\t"+breaks.get(breaks.size()-1)+"\t"+exonCount+
				"\t"+CigarHash2.getString(breaks)+"\t"+exon_str+"\t"+cnt.sum()+"\t"+read_count1+"\t";
				o.printTranscript(str,depth_str);
			}
		//	o.printTranscriptAlt(cc);
			//	String geneP_header = "ID\tchrom\tstart\tend\tgene_nme\tnum_exons\tisoforms\tbreaks\ttotLen\tcountTotal\t"+TranscriptUtils.getString("count", num_sources,true);

			String str1 = cc.id()+"\t"+chrom+"\t"+cc.startPos+"\t"+cc.endPos+"\t"+
					cc.exonCount()+"\t"+cc.numIsoforms()+"\t"+CigarHash2.getString(new ArrayList<String>(genes))+"\t"+
					"\t"+cc.readCountSum()+"\t"+read_count;
			o.printGene(str1,depth_str);
			List<Integer> br = cc.getBreaks();
		
			StringBuffer starts = new StringBuffer();
			StringBuffer ends = new StringBuffer();
			StringBuffer chroms = new StringBuffer();
			StringBuffer strand = new StringBuffer();
				int len =0;
				if(br!=null){
				for(int i=0; i<br.size(); i+=2){
					len += br.get(i+1) - br.get(i);
					String ext = i<br.size()-2 ? ";" : "";
					starts.append(br.get(i)+ext); ends.append(br.get(i+1)+ext);
					strand.append("+"+ext); chroms.append(chrom+ext);
				}
				}else{
					starts.append(cc.startPos); ends.append(cc.endPos); chroms.append(cc.endPos); strand.append("+");
				}
				
			
			o.printFC(cc.breaks_hash.secondKey+"\t"+chroms.toString()+"\t"+starts.toString()+"\t"+ends.toString()+"\t"+"+;+\t"+len+"\t"+read_count);
		//	Geneid  Chr     Start   End     Strand  Length  /DataOnline/Data/Projects/corona_invitro/host_analysis/direct_cDNA/vero/vero_24hpi/merged/genome/infected/mo
		//	ID0.0   MT007544.1;MT007544.1   14;27385        66;29860        +;+     2529    0       0       0       0       0       0
			//I
//					CigarCluster.recordDepthByPosition ?  cc.getTotDepthSt(true)+"\t"+cc.getTotDepthSt(false): "");
		}

		public String getError(int source_index) {
			if(this.map==null) return "";
			else return this.map.getError(source_index);
		}

		

		

		


		
		
		
	}