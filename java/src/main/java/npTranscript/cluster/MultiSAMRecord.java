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

import htsjdk.samtools.AlignmentBlock;
import htsjdk.samtools.SAMRecord;
import npTranscript.NW.NeedlemanWunsch;
import npTranscript.run.TranscriptUtils1;

public class MultiSAMRecord{
	
	static class StringBuffers{
		StringBuffer chrom = new StringBuffer();
		StringBuffer strand = new StringBuffer();//sam1.get(0).getReadNegativeStrandFlag() ? '-': '+';
		StringBuffer q_str = new StringBuffer();
		StringBuffer q_str1 = new StringBuffer();
		StringBuffer read_pos = new StringBuffer();
		StringBuffer ref_pos = new StringBuffer();
		StringBuffer ref_3 = new StringBuffer();
		StringBuffer overlaps = new StringBuffer();
		public void append(String string) {
			read_pos.append(";");
			chrom.append(";");//strand.append(";");
			q_str.append(";");q_str1.append(";");
			overlaps.append(";");
			ref_pos.append(";");
			ref_3.append(";");
			// TODO Auto-generated method stub
			
		}
		public String getString() {
			// TODO Auto-generated method stub
			return chrom.toString()+" "+strand.toString()+" "+ref_3.toString()+" "+ q_str.toString()+" "+q_str1.toString()+" " +read_pos.toString()+" "+ref_pos.toString()+" "+overlaps.toString();
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


	public static String toString1(List<Integer> pos) {
		StringBuffer sb = new StringBuffer();
		for(int j=0; j<pos.size(); j+=2) {
			if(j>0) sb.append(",");
			sb.append(pos.get(j)+"_"+pos.get(j+1));
		}
		return sb.toString();
	}
public static class Breaks{
	List<Integer> pos = new ArrayList<Integer>();
	boolean neg1;
	boolean prime3;
	Breaks(SAMRecord sr, int gap, boolean prime3){
		this.prime3= prime3;
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
		if(!neg1 && prime3 || neg1 && !prime3) {
			this.pos = this.pos.reversed();
		}
	}
	 public String first() {
		 return pos.get(0).toString();
	 }
	public String toString() {
		return(toString1(this.pos));
	}
	
}
	
		//SAMRecord s1;
		String read_str, readname;
		int read_length;
		boolean primary_neg, prime3; 
	//	Integer[] read_coords = Integer[];
		byte[] bq;
		
		StringBuffers sb = new StringBuffers();
		// len; 
		TreeMap<Integer[], SAMRecord> res;
		TreeMap<Integer[], Breaks> res1;
		int gap;
		
		MultiSAMRecord(SAMRecord primary, int gap, boolean prime3){
			this.gap = gap;
			this.prime3=prime3;
			this.readname = primary.getReadName();
			bq = primary.getBaseQualities();
			if(primary.isSecondaryAlignment()) throw new RuntimeException("!");
		//	this.s1 = primary;
			read_str = primary.getReadString();
			read_length = primary.getReadLength();
			primary_neg = primary.getReadNegativeStrandFlag();
			//this.primary=true;
			res = new TreeMap<Integer[],SAMRecord>(prime3 ? reverse_comp : forward_comp);
			res1 = new TreeMap<Integer[],Breaks>(prime3 ? reverse_comp : forward_comp);
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
				Breaks br = new Breaks(supp, gap, prime3);
				res1.put(key, br);
			}
		}
		
		public String toString() {
			StringBuffer sb = new StringBuffer(this.readname+"\n");
			sb.append(IdentityProfileHolder.gson.toJson(res.keySet()));
			sb.append(IdentityProfileHolder.gson.toJson(copy(res).values()));
			
			return(sb.toString());
		}

		
		
		
		public boolean checkList(	){
			Integer[] se_prev=null;
			boolean removed=false;
			String read_str1 = primary_neg ?  NeedlemanWunsch.reverseComplement(read_str) : read_str;
			int cnt=1;
			Set<Integer[]> keyset = prime3 ? res.descendingKeySet() : res.keySet();
			Iterator<Integer[]> it1 = keyset.iterator();
			
			
			
			while(it1.hasNext()) {
					boolean first = se_prev==null;
				
					Integer[] se1 = it1.next();
					boolean last = !it1.hasNext();
					SAMRecord sr2 = res.get(se1);
					String char2 = sr2.getReadNegativeStrandFlag() ? "-" : "+";
					int st2 = sr2.getStart();
					int len = sr2.getReadLength();
					if(!first) sb.append(";");
					sb.strand.append(sr2.getReadNegativeStrandFlag() ? '-': '+');
					
					sb.q_str.append(sr2.getMappingQuality());
					//Breaks br2 = res1.get(se1);
					
					List<Integer> read_pos = this.getReadPos(se1);
					List<Integer> ref_pos = this.getRefPos(se1);
					boolean neg = sr2.getReadNegativeStrandFlag();
					String header = ">"+readname+" cnt:"+cnt+" "+char2+" "+sr2.getReferenceName()+" "+st2;
					if(Outputs.spliceOut!=null && read_pos.size()>2) {
						for(int j=2; j<read_pos.size(); j+=2) {
							int pos1 = read_pos.get(j-1);
							int pos2 = read_pos.get(j);
						
							String substr = (neg && ! prime3 || !neg && prime3)  ? read_str1.substring(pos2-3, pos1+2) : read_str1.substring(pos1-3, pos2+2);
							String substr1 = (neg && ! prime3 || !neg && prime3)  ? read_str1.substring(pos2-2, pos1+1) : read_str1.substring(pos1-2, pos2+1);
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
					
					
					sb.ref_pos.append(toString1(ref_pos));
					sb.ref_3.append(ref_pos.getFirst().toString());
					
					int qstart = se1[0];//record.getReadPositionAtReferencePosition(record.getAlignmentStart());
					int qend = se1[1];//record.getReadPositionAtReferencePosition(record.getAlignmentEnd());
					sb.read_pos.append(toString1(read_pos));
					sb.q_str1.append(TranscriptUtils1.median(bq, qstart, qend-qstart));
				 	sb.chrom.append(sr2.getReferenceName());
					
					
					if(!first) {
						sb.append(";");
						int overlap = se_prev[1] - se1[0]+1;
						sb.overlaps.append(overlap);
						SAMRecord sr1 = res.get(se_prev);
						//SAMRecord sr2 = res.get(se1);
						String char1 = sr1.getReadNegativeStrandFlag() ? "-" : "+";
						int st1 = sr1.getStart();
						boolean same =sr1.getReferenceName().equals(sr2.getReferenceName());
						String header1 = ">"+readname+",overlap:"+overlap+" cnt:"+cnt+" "+se1[0]+","+se_prev[1]+" "+char1+"," +char2+" "+sr1.getReferenceName()+","+sr2.getReferenceName()+" "+st1+" "+st2+" "+same;
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
								String substr = read_str1.substring(se_prev[1]-1, se1[0]-1);
								//if(primary_neg) substr =  NeedlemanWunsch.reverseComplement(substr);
								Outputs.joinOut.println(header1+"\t"+substr);
							}
							if(overlap < -1 *IdentityProfileHolder.gap_max) {
								removed=true;
							}
						}else if(overlap<2 && overlap >-2) {
							if(Outputs.noGap!=null) {
								String substr = read_str1.substring(se_prev[1]-6, se1[0]+5); //because its one based
								
								String substr1 = read_str1.substring(se_prev[1]-4, se1[0]+3); //because its one based
								String substr2 = read_str1.substring(se_prev[1]-3, se1[0]+2); 
								//if(primary_neg) substr =  NeedlemanWunsch.reverseComplement(substr);
								Outputs.noGap.println(header1+"\t"+substr+"\t"+substr1+"\t"+substr2);
							}
							
						}
						cnt++;
					}
					
					
					
					
					se_prev = se1;
					
					
					
				}
				
				if(Outputs.allOut!=null) {
					Outputs.allOut.println(">"+readname+" "+sb.getString()+"\n"+read_str1);
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