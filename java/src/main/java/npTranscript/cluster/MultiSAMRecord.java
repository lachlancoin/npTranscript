package npTranscript.cluster;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.Iterator;
import java.util.List;
import java.util.Map.Entry;
import java.util.Set;
import java.util.SortedMap;
import java.util.TreeMap;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import htsjdk.samtools.AlignmentBlock;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMTag;
import npTranscript.NW.NeedlemanWunsch;
import npTranscript.run.TranscriptUtils1;

public class MultiSAMRecord{
	
	public static String toString1(List<Integer> pos, double round) {
		StringBuffer sb1 = new StringBuffer();
		for(int j=0; j<pos.size(); j+=2) {
			if(j>0) sb1.append(",");
			sb1.append(Math.round( pos.get(j).doubleValue()/round)+"_"+Math.round(pos.get(j+1).doubleValue()/round));
		}
		return sb1.toString();
	}
	public static String toString2(List<List<Integer>> pos, double round) {
		StringBuffer sb1 = new StringBuffer();
		for(int j=0; j<pos.size(); j++) {
			if(j>0) sb1.append(";");
			sb1.append(toString1(pos.get(j), round));
		}
		return sb1.toString();
	}
	
	public static String toString3(List<String> l , boolean flip, CharSequence sep) {
		List<String> l1 = flip ? l.reversed() : l;
		return l1.stream().collect(Collectors.joining(sep));
	}
	
	static class StringBuffers{
		String readname;
		boolean prime3=false;
		Integer polyA;
		StringBuffers(){
			
		}
		void reset(String readname, Integer polyA){
			this.readname = readname;
			this.polyA = polyA;
			this.prime3 = false;
		}
		public void set3prime(boolean prime3) {
			this.prime3 = prime3;
		}
		List<String >chrom = new ArrayList<String>();
		List<String>strand = new ArrayList<String>();//sam1.get(0).getReadNegativeStrandFlag() ? '-': '+';
		List<String> q_str = new ArrayList<String>();
		List<String> q_str1 = new ArrayList<String>();
		List<List<Integer>> ref_pos = new ArrayList<List<Integer>>();
		List<List<Integer>> read_pos = new ArrayList<List<Integer>>();
		//StringBuffer ref_3 = new StringBuffer("e:");
		List<String> overlaps = new ArrayList<String>();
		
		public String getString() {
			return getString(1);
			
		}
		public String chrom() {return toString3(chrom, prime3, ";");}
		public String strand() {
			return toString3(strand, prime3, ";");
			}
		public String q_str() {return toString3(q_str, prime3, ";");}
		public String q_str1() {return toString3(q_str1, prime3, ";");}
		public String overlaps() {return toString3(overlaps, prime3, ";");}
		
		public String end(double round) {
			String res = this.ref_pos.stream().map(t->String.valueOf((int) Math.round(t.get(0).doubleValue()/round))).collect(Collectors.joining(";"));
			return res;
		}
		public String ref(double round) {
			return toString2(this.ref_pos, round);
		}
		public String read(double round) {
			return toString2(this.read_pos, round);
		}
		
		
		public String getString(double round) {
			// TODO Auto-generated method stub
			return chrom()+" "+strand()+" "+end(round)+" "+ q_str()+" "+q_str1()+" " +read(round)+" "+ref(round)+" "+overlaps();
		}
		
		public void append(String chrom, char strand, int q_str, double q_str1, List<Integer> ref_pos2,List<Integer> read_pos2) {
			this.chrom.add(chrom);
			this.strand.add(strand+"");
			this.q_str.add(q_str+"");
			this.q_str1.add(q_str1+"");
		//	this.overlaps.add(overlaps);
			this.ref_pos.add(ref_pos2);
			this.read_pos.add(read_pos2);
			// TODO Auto-generated method stub
			
		}
		public Integer polyA() {
			// TODO Auto-generated method stub
			return polyA;
		}
	}
	
	
	
	static class ReadOrderComparator implements Comparator<Integer[]>{
		boolean prime3;
		
		ReadOrderComparator(boolean prime3){
			this.prime3 = prime3;
		}
		public void setPrime3(boolean val) {
			this.prime3 = val;
		}
		
		@Override
		public int compare(Integer[] se1, Integer[] se2) {
			if(se1[0]==se2[0]) return 0;
			if(prime3) {
				return se1[0] < se2[0]? 1 : -1;
			}else {
				return se1[0] < se2[0]? -1 : 1;
			}
		}

		
	}
	
	static ReadOrderComparator forward_comp = new ReadOrderComparator(false);
	static ReadOrderComparator reverse_comp = new ReadOrderComparator(false);
	
	
	static SortedMap<Integer[], String> copy(SortedMap<Integer[], SAMRecord> map) {
		TreeMap<Integer[], String> tm = new TreeMap<Integer[], String>(map.comparator());
		Iterator<Entry<Integer[], SAMRecord>> it = map.entrySet().iterator();
		while(it.hasNext()) {
			Entry<Integer[], SAMRecord> ent = it.next();
			SAMRecord sr = ent.getValue();
			Integer[] key = ent.getKey();
			tm.put(key,key[0]+","+key[1]+":"+
					sr.getReferenceName()+",neg:"+sr.getReadNegativeStrandFlag()+",supp:"+sr.isSecondaryOrSupplementary());
		}
		return tm;
	}


	
public static class Breaks{
	List<Integer> pos = new ArrayList<Integer>();
	boolean neg1;
	Breaks(SAMRecord sr, int gap){
		Iterator<AlignmentBlock> ab = sr.getAlignmentBlocks().iterator();
		AlignmentBlock ab1 = ab.next();
			int st = ab1.getReferenceStart();
			int end = st + ab1.getLength()-1;
			
			while(ab.hasNext()) {
				ab1=ab.next();
				int st1 = ab1.getReferenceStart();
				if(st1 - end > gap) {
					pos.add(st); pos.add(end);
			//		int a1 =sr.getReadPositionAtReferencePosition(st);
				//	int a12 =sr.getReadPositionAtReferencePosition(end);
					st = st1;
				}
					end = st1 +ab1.getLength()-1;
			}
			if(pos.contains(st)) throw new RuntimeException("!!");
			pos.add(st); pos.add(end);
			
		 neg1 = sr.getReadNegativeStrandFlag();
		//if neg1 is false then it will be in 5 -> 3
		if(neg1) {
			this.pos = this.pos.reversed();
		}
	}
	 public String first() {
		 return pos.get(0).toString();
	 }
	public String toString() {
		return(toString1(this.pos,1.0));
	}
	
}
	
		//SAMRecord s1;
		String read_str, readname;
		int read_length;
		boolean primary_neg;//, prime3; 
	//	Integer[] read_coords = Integer[];
		byte[] bq;
		
		StringBuffers sb;
		// len; 
		TreeMap<Integer[], SAMRecord> res;
		TreeMap<Integer[], Breaks> res1;
		int gap;
		
		MultiSAMRecord(int gap){
			this.gap = gap;
			this.sb  = new StringBuffers();
			res = new TreeMap<Integer[],SAMRecord>(forward_comp);
			res1 = new TreeMap<Integer[],Breaks>(forward_comp);
		}
		
		public void reset(SAMRecord primary){
			Integer polyA = (Integer) primary.getAttribute(SAMTag.PT);
			this.readname = primary.getReadName();
			sb.reset(readname, polyA);
			bq = primary.getBaseQualities();
			if(primary.isSecondaryAlignment()) throw new RuntimeException("!");
		//	this.s1 = primary;
			read_str = primary.getReadString();
			read_length = primary.getReadLength();
			primary_neg = primary.getReadNegativeStrandFlag();
			//this.primary=true;
			res1.clear();res.clear();
			this.update(primary, true);
		}
		
		public void update(SAMRecord supp, boolean add) {
			//boolean primary = !supp.isSecondaryOrSupplementary();
			int len = supp.getReadLength();
			String str1 = supp.getReadString();
			boolean neg1=supp.getReadNegativeStrandFlag();
		//	this.s1 = supp;
			if(add) {
				Integer[] key =this.getStEnd(supp); 
				res.put(key, supp);
				Breaks br = new Breaks(supp, gap);
				res1.put(key, br);
			}
		}
		
		public String toString() {
			StringBuffer sb1 = new StringBuffer(this.readname+"\n");
			sb1.append(IdentityProfileHolder.gson.toJson(res.keySet()));
			sb1.append(IdentityProfileHolder.gson.toJson(copy(res).values()));
			
			return(sb1.toString());
		}

		
		
		
		public boolean checkList(	){
			Integer[] se_prev=null;
			boolean removed=false;
			String read_str1 = primary_neg ?  NeedlemanWunsch.reverseComplement(read_str) : read_str;
			int cnt=1;
			Set<Integer[]> keyset =  res.keySet();//prime3 ? res.descendingKeySet() :
			Iterator<Integer[]> it1 = keyset.iterator();
			
			
			
			while(it1.hasNext()) {
					boolean first = se_prev==null;
				
					Integer[] se1 = it1.next();
					boolean last = !it1.hasNext();
					SAMRecord sr2 = res.get(se1);
					String char2 = sr2.getReadNegativeStrandFlag() ? "-" : "+";
					int st2 = sr2.getStart();
					int len = sr2.getReadLength();
					
					//Breaks br2 = res1.get(se1);
					
					List<Integer> read_pos = this.getReadPos(se1);
					List<Integer> ref_pos = this.getRefPos(se1);
					boolean neg = sr2.getReadNegativeStrandFlag();
					String header = ">"+readname+" cnt:"+cnt+" S"+char2+" "+sr2.getReferenceName()+" "+st2;
					if(Outputs.spliceOut!=null && read_pos.size()>2) {
						for(int j=2; j<read_pos.size(); j+=2) {
							int pos1 = read_pos.get(j-1);
							int pos2 = read_pos.get(j);
						
							String substr = (neg)  ? read_str1.substring(pos2-3, pos1+2) : read_str1.substring(pos1-3, pos2+2);
							String substr1 = (neg )  ? read_str1.substring(pos2-2, pos1+1) : read_str1.substring(pos1-2, pos2+1);
								Outputs.spliceOut.println(header+" "+pos1+"-"+pos2+"\t"+substr+"\t"+substr1);
						}
						
					}
					int firstPos = read_pos.get(0);
					int lastPos = read_pos.get(read_pos.size()-1);
					if(first && Outputs.fiveOut!=null && (!neg && firstPos==1 || neg && firstPos==len)) {
						String substr=neg ? read_str.substring(firstPos-10, firstPos) :  read_str.substring(firstPos-1, firstPos+9);
						String substr1=neg ?read_str.substring(firstPos-5, firstPos) : read_str.substring(firstPos-1, firstPos+4);
						if(neg) {
							substr = NeedlemanWunsch.reverseComplement(substr);
							substr1 = NeedlemanWunsch.reverseComplement(substr1);
						}
						
//						String substr =  read_str1.substring(0, 10);
//						String substr1 = read_str1.substring(0, 5);
						Outputs.fiveOut.println(header+"\t"+substr+"\t"+substr1);
					}
					if(last && Outputs.threeOut!=null && (!neg && lastPos>len-100 || neg && lastPos<100)) {
						String substr=neg ? read_str.substring(lastPos-1, lastPos+9): read_str.substring(lastPos-10, lastPos);
						String substr1=neg ? read_str.substring(lastPos-1, lastPos+4): read_str.substring(lastPos-5, lastPos);
						if(neg) {
							substr = NeedlemanWunsch.reverseComplement(substr);
							substr1 = NeedlemanWunsch.reverseComplement(substr1);
						}
						
						//String substr =  read_str1.substring(len-10, len);
						//String substr1 = read_str1.substring(len-5, len);
						Outputs.threeOut.println(header+"\t"+substr+"\t"+substr1);
					}
					
					//sb.ref_pos.append(toString1(ref_pos));
					
					int qstart = se1[0];//record.getReadPositionAtReferencePosition(record.getAlignmentStart());
					int qend = se1[1];//record.getReadPositionAtReferencePosition(record.getAlignmentEnd());
					//sb.appendRead(read_pos);
				//	sb.read_pos.append(toString1(read_pos));
					
					sb.append(sr2.getReferenceName(),sr2.getReadNegativeStrandFlag() ? '-': '+',sr2.getMappingQuality(),
							TranscriptUtils1.median(bq, qstart, qend-qstart), ref_pos,read_pos);
				//	public void append(String chrom, String strand, String q_str, String q_str1, String overlaps, List<Integer> ref_pos2,List<Integer> read_pos2) {

					if(!first) {
						int overlap = se_prev[1] - se1[0]+1;
						sb.overlaps.add(overlap+"");
						SAMRecord sr1 = res.get(se_prev);
						//SAMRecord sr2 = res.get(se1);
						String char1 = sr1.getReadNegativeStrandFlag() ? "-" : "+";
						int st1 = sr1.getStart();
						boolean same =sr1.getReferenceName().equals(sr2.getReferenceName());
						String header1 = ">"+readname+",overlap:"+overlap+" cnt:"+cnt+" "+se1[0]+","+se_prev[1]+" S"+char1+char2+" "+sr1.getReferenceName()+","+sr2.getReferenceName()+" "+st1+" "+st2+" "+same;
						if(overlap>0) {
							if(Outputs.overlapOut!=null) {
								String substr = read_str1.substring(se1[0]-1,se_prev[1]-1);
							//	if(primary_neg) substr =  NeedlemanWunsch.reverseComplement(substr);
								Outputs.overlapOut.println(header1+"\t"+substr);
							}
							if(overlap>IdentityProfileHolder.overlap_max) {
								removed=true;
							}
							//return null;

						}else if(overlap <0) {
							if(Outputs.joinOut!=null) {
								String substr = read_str1.substring(se_prev[1]-1, se1[0]);
								//if(primary_neg) substr =  NeedlemanWunsch.reverseComplement(substr);
								Outputs.joinOut.println(header1+"\t"+substr);
							}
							if(overlap < -1 *IdentityProfileHolder.gap_max) {
								removed=true;
							}
						}else if(overlap<2 && overlap >-2) {
							if(Outputs.noGap!=null) {
								String substr1 = read_str1.substring(se_prev[1]-6, se_prev[1]); //because its one based
								String g = overlap>0 ? read_str1.substring(se_prev[1], se1[0]) : "";
								String substr2 = read_str1.substring(se1[0]-1, se1[0]+6); //because its one based
								
//								String substr1 = read_str1.substring(se_prev[1]-4, se1[0]+3); //because its one based
//								String substr2 = read_str1.substring(se_prev[1]-3, se1[0]+2); 
								//if(primary_neg) substr =  NeedlemanWunsch.reverseComplement(substr);
								Outputs.noGap.println(header1+"\t"+substr1+"\t"+g+"\t"+substr2);
							}
							
						}
						cnt++;
					}
					
					
					
					
					se_prev = se1;
					
					
					
				}
				
				if(Outputs.allOut!=null) {
					String head = sb.getString();
					Outputs.allOut.println(">"+readname+" "+head+"\t"+read_str1);
				}
				
				
				return removed;
				
		}
		
		public List<Integer> getReadPos(Integer[] key){
			SAMRecord sr = this.res.get(key);
			List<Integer> ref_br_pos = this.res1.get(key).pos;
			List<Integer>readpos = (toRead(sr,ref_br_pos));
			return(readpos);
		}
		public List<Integer> getRefPos(Integer[] key){
			List<Integer>readpos = ((this.res1.get(key).pos));
			return(readpos);
		}
		
		public Integer[] getStEnd(SAMRecord s1) {
			Integer[] ref = new Integer[] {s1.getStart(), s1.getEnd()};
			
			List<Integer> reads = toRead(s1, Arrays.asList(ref));
			return(new Integer[] {reads.get(0),reads.get(1)});
		}
		public List<Integer> toRead(SAMRecord s1, List<Integer> ref) {
		//	Integer[] read_coords = new int[2];
			String read_str1 = this.read_str;
            String str1 = s1.getReadString();
            boolean neg1 = s1.getReadNegativeStrandFlag();
            int len = s1.getReadLength();
			List<Integer> reads_ = ref.stream().map(t-> s1.getReadPositionAtReferencePosition(t)).collect(Collectors.toList());

           
//			if(!neg1 & st <end & st_>end_) {
//				throw new RuntimeException("this should not happen");
	//		}
			boolean adjusted=false;
			if(len!=read_length   ) {
		//		if(st_>5) throw new RuntimeException("!! should not happen");
				String str = str1;
				if(neg1) { 
					// System.err.println(str);
					str = NeedlemanWunsch.reverseComplement(str1);
					// System.err.println(str);
				}
				if(primary_neg) {
					read_str1 = NeedlemanWunsch.reverseComplement(read_str);
				}
				
				int shift_read = read_str1.indexOf(str);
				if(shift_read<0) {
					throw new RuntimeException("no subindex!");
				}
				//  st_= shift_read+1;
				 // end_ = st_ + len;
				if(shift_read>0) {
				reads_ = reads_.stream().map(t-> t+shift_read).collect(Collectors.toList());
				}
				 adjusted=true;
			}
			  if(neg1 & !adjusted) { // & !primary_neg || primary_neg && !neg1) {
			    reads_ = reads_.stream().map(t -> read_length -t+1).collect(Collectors.toList()).reversed()	;
			    }
			//  System.err.println(reads_);
			  return(reads_);
		}

		public boolean fusion() {
			return(this.res.size()>1);
		}
		
	}