package npTranscript.cluster;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.Iterator;
import java.util.List;
import java.util.Map.Entry;
import java.util.SortedMap;
import java.util.SortedSet;
import java.util.Stack;
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.ThreadPoolExecutor;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.stream.Collectors;

import com.google.gson.Gson;

import htsjdk.samtools.AlignmentBlock;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.SequenceUtil;
import japsa.seq.Alphabet;
import japsa.seq.Sequence;
import npTranscript.NW.AlignmentParameters;
import npTranscript.NW.AlignmentResult;
import npTranscript.NW.NeedlemanWunsch;
import npTranscript.NW.PolyAT;
import npTranscript.run.TranscriptUtils1;
import npTranscript.run.ViralTranscriptAnalysisCmd2;

public class IdentityProfileHolder {
	
	static boolean CHECK=true;
	public static double overlap_thresh = 0.33;
	public static double overlap_max = 20;
	public static double gap_max = 20;  // max gap on read between adjacent alignments

	static Gson gson = new Gson();
	public static void waitOnThreads(int sleep) {
		if(executor==null) return;
		if(executor instanceof ThreadPoolExecutor){
	    	while(((ThreadPoolExecutor)executor).getActiveCount()>0){
	    		try{
		    	//System.err.println("IdentityProfileHolder: awaiting completion "+((ThreadPoolExecutor)executor).getActiveCount());
		    	//Thread.currentThread();
				Thread.sleep(sleep);
	    		}catch(InterruptedException exc){
	    			exc.printStackTrace();
	    		}
	    	}
	    	}
		//executor.a
	}
	
	static int primelen = 500;//chr.length();
	public static  ExecutorService executor ;
 //final CigarClusters all_clusters;
  BreakPoints bp;
 final Outputs o;
 //Sequence genome;
 // String chrom_="";
 // int chrom_index;
 final String type_nmes;
 final int num_sources;
 //Sequence chr5prime; Sequence chr3prime;
 
 Stack<IdentityProfile1> idents = new Stack<IdentityProfile1>();
 
 
 public synchronized void replace(IdentityProfile1 pr1){
	 //pr1.clear();
	// System.err.print("returning");
	 idents.push(pr1);
	 
 }
 /** supp is in case this is a supplementary alignment, in which case we need to get the same profile as previously used 
  * we assume that the primary alignment is first, and supplementary alignments follow
  * note that this breaks multi-threading
  * */
 public synchronized IdentityProfile1 get(){
	 int len = idents.size();
	 IdentityProfile1 idp;
	 if(len==0){
		idp =  new IdentityProfile1(this) ;
	 }else{
		 idp = idents.pop();
	
		 idp.clear(); // gets object clear for next use

	 }
	return idp;
 }
 
 public static int cleanUpMoreThan = 100000; //  distance between previous end and current start to remove cluster
	public static int clearUpThreshold = 100; //number of clusters before try to clear up 
	
	/*
 public void clearUpTo(int endThresh){
	 
	 this.all_clusters.clearUpTo(endThresh, o, genome, chrom_index);
 }
 public void printBreakPoints(int chrom_index) throws IOException {
		if(bp!=null) this.bp.printBreakPoints(this.o, chrom_index);
		
	}
*/ 
 
	/*public void updateChrom(Sequence refSeq, String chrname, int chrom_index) {
		// TODO Auto-generated method stub
		if(this.chrom_.equals(chrname)) return;
		 this.genome=refSeq;
		 this.chrom_ =chrname;
		 this.chrom_index = chrom_index;
		 chr5prime = TranscriptUtils.coronavirus ? refSeq.subSequence(0	, Math.min( primelen, refSeq.length())) : null;
		chr3prime = TranscriptUtils.coronavirus ? refSeq.subSequence(Math.max(0, refSeq.length()-primelen), refSeq.length()) : null;
		if(calcBreakpoints && genomes!=null){
			System.err.println("calculating break point usage");
			bp = new BreakPoints(num_sources, refSeq, genomes);
		}else{
			bp  = null;
		}
	}*/
 
	boolean calcBreakpoints;
	ArrayList<Sequence> genomes;
	Annotation annot;
 public IdentityProfileHolder(ArrayList<Sequence> genomes, Sequence refSeq,
		 Outputs o, String in_nmes,  boolean calcBreakpoints, Annotation annot)
		 throws IOException {
	 this.genomes = genomes;
	 this.annot = annot;
	 this.calcBreakpoints = calcBreakpoints;
	 Arrays.fill(basesStart,0);
	 Arrays.fill(basesEnd,0);
	 this.type_nmes = in_nmes;
	 this.num_sources = 1;//in_nmes.length;
	 this.o  = o;
	// all_clusters =new CigarClusters( num_sources, genomes);
	 bp=null;
	 if(calcBreakpoints && genomes!=null){
			System.err.println("calculating break point usage");
			bp = new BreakPoints(num_sources, refSeq, genomes);
	 }
	}
 
 static int lengthOnRead(SAMRecord s1) {
		List<AlignmentBlock> l1 = s1.getAlignmentBlocks();
		AlignmentBlock a1 = l1.get(0);
		AlignmentBlock a2 = l1.get(l1.size()-1);
		int st1 = a1.getReadStart();
		int end1 = a2.getReadStart()+a2.getLength();
		int diff = end1 -st1;
		return diff;
	}
 static int readStart(SAMRecord s1) {
		List<AlignmentBlock> l1 = s1.getAlignmentBlocks();
		AlignmentBlock a1 = l1.get(0);
		AlignmentBlock a2 = l1.get(l1.size()-1);
		int st1 = a1.getReadStart();
		int end1 = a2.getReadStart()+a2.getLength();
		int diff = end1 -st1;
		return st1;
	}
static class SamComparator implements Comparator<SAMRecord>{
	boolean quality;
	
	public SamComparator(boolean b) {
		this.quality = b;
	}

	@Override
	public int compare(SAMRecord o1, SAMRecord o2) {
		if(quality) {
			int i1 = Boolean.compare(o1.isSecondaryOrSupplementary(), o2.isSecondaryOrSupplementary());
			if(i1!=0) return i1;
			int i2 = -1*Integer.compare(o1.getMappingQuality(), o2.getMappingQuality());
			if(i2!=0) return i2;
			
			
			return 1*Integer.compare(readStart(o1),readStart(o2));
		}
		else return  Integer.compare(o1.getAlignmentBlocks().get(0).getReadStart(),
			o2.getAlignmentBlocks().get(0).getReadStart());
	}
	
};
static Comparator comp = new SamComparator(false);
static Comparator comp_q = new SamComparator(true);

static AlignmentParameters align_p =new AlignmentParameters();
	
	static SAMRecord groupSam(List<SAMRecord>sam_1, List<SAMRecord> extracted,
			List<int[]>se){
		Iterator<SAMRecord> sams = sam_1.iterator();
		SAMRecord first = sams.next();
		int mq = first.getMappingQuality();
		SAMRecord primary  = null;
		if(!first.isSecondaryOrSupplementary()) primary = first;
		if(primary==null) throw new RuntimeException("first is not primary");
		int read_length = primary.getReadLength();
		String read_str = primary.getReadString();
		int q1 = first.getMappingQuality();
		//sams.remove();
		extracted.add(first);
		boolean primary_neg = primary.getReadNegativeStrandFlag();
	
		//int[] overl = new int[3];
		if(sams.hasNext()) {
			se.add(getStEnd(first, read_str,read_length, primary_neg));
			outer: for(int i=0; sams.hasNext(); i++){
				SAMRecord sam = sams.next();
				//int mq1 = sam.getMappingQuality();
				int[] se1=getStEnd(sam, read_str,read_length, primary_neg);
				for(int j=0; j<extracted.size(); j++){
					int[] se_j=se.get(j);
					int overl1 = read_overlap(se1,se_j);
					System.err.println(overl1);
					if( overl1>overlap_thresh) {
					     continue outer;
					}
				}
					extracted.add(sam);
					se.add(se1);
			}
		}
		
		//Collections.sort(extracted , comp);
		
		/*if(false){
		for(int i=0; i<extracted.size(); i++){
			SAMRecord sr = extracted.get(i);
			System.err.println(sr.getReferenceName() + " "+sr.getReadNegativeStrandFlag());
			AlignmentBlock ab1 = sr.getAlignmentBlocks().get(0);
			AlignmentBlock ab2 = sr.getAlignmentBlocks().get(sr.getAlignmentBlocks().size()-1);
			System.err.println(ab1.getReadStart()+" - "+ab1.getReferenceStart());
			System.err.println(ab2.getReadStart()+" - "+ab2.getReferenceStart());

		}
		System.err.println("h");
		}*/
	//	return extracted;
		//return extracted;
		//if(flipped!=null) return flipped.intValue()==0 ? fa
		return primary;
	}
	
	static Pattern patt = Pattern.compile("[ACTG]");
 
	static int[] getStEnd(SAMRecord s1,  String read_str, int read_length, boolean primary_neg) {
		boolean primary= !s1.isSecondaryOrSupplementary();
		int mq = s1.getMappingQuality();
		String str1 =s1.getReadString(); String read_str1 = read_str;
		int len =s1.getReadLength();
		boolean neg1=s1.getReadNegativeStrandFlag();

		
		//int len1 = s1.getAlignmentEnd()-s1.getAlignmentStart();
		//List<AlignmentBlock> l1 = s1.getAlignmentBlocks();
		//AlignmentBlock a1 = l1.get(0);
	//	AlignmentBlock a2 = l1.get(l1.size()-1);
		int st = s1.getStart();int end = s1.getEnd();
		int st_ = s1.getReadPositionAtReferencePosition(st);
		int end_ = s1.getReadPositionAtReferencePosition(end);
		if(!neg1 & st <end & st_>end_) {
			throw new RuntimeException("this should not happen");
		}
//	    if(neg1 && ! primary_neg || !neg1 && primary_neg) {
		boolean aligns_to_end=false;
		boolean adjusted=false;
		if(len!=read_length   & st_<5 & true) {
			String str = str1;
			if(neg1) { 
				// System.err.println(str);
				str = NeedlemanWunsch.reverseComplement(str1);
				// System.err.println(str);
			}
			if(primary_neg) {
				read_str1 = NeedlemanWunsch.reverseComplement(read_str);
			}
			
			int shift_read = read_str1.indexOf(str)+1;
			if(shift_read<0) {
				if(true) throw new RuntimeException("no subindex!");
					AlignmentResult result= 
							NeedlemanWunsch.computeNWAlignment(read_str1, str, align_p);
					double sc = NeedlemanWunsch.score(result);
					  String[] alignments = result.getAlignments();
					
					//  System.err.println("ref\n"+alignments[0]);
					  String s2 = alignments[1];
					  Matcher m = patt.matcher(s2);
					  m.find();
					   shift_read = m.start();
			}
			 
		//	  st_ =  st_+shift_read;
			//  end_ = end_+shift_read;
			  st_ = shift_read;
			  end_ = st_ + len;
			 adjusted=true;
		}
		
		
	    if(neg1 & !adjusted) { // & !primary_neg || primary_neg && !neg1) {
	    	
	    	int st1 = read_length-end_ + 1;
	    	int st2 = read_length - st_ +1;
	    //	System.err.println(read_str.substring(st1,st2));
	    //	System.err.println(read_str.substring(st_,end_));
	    //	System.err.println(NeedlemanWunsch.reverseComplement(str1));
	    	return new int[] {st1,st2};
	    }
		
		return new int[] {st_,end_};
/*		int st1 = a1.getReadStart();
		int end1 = a2.getReadStart()+a2.getLength();
	//	int readlen = s1.getReadLength();
		if(neg1) {
			int end2 = read_length-st1+1;
			int st2 = read_length-end1+1;
			
			 
			
//System.err.println(s1+" "+st2+" "+end2);
			return new int[] {st2,end2};

		}
//		System.err.println(s1+" "+st1+" "+end1);

		return new int[] {st1,end1};*/
	}

	private static int read_overlap(int[] se1, int[] se2) {
		int st1 = se1[0]; int end1 = se1[1];
		int st2 = se2[0]; int end2 = se2[1];
	
		int overlap = CigarHash2.overlap(st2, end2,st1 , end1);
	/*	if(overlap>overl[2] || overl[2] ==0){
			overl[0] = end2 -st2; 
			overl[1] = end1 -st1;
			overl[2] = overlap;
		}*/
		return overlap;
	}
 
	
	public synchronized void addBreakPoint(int source_index, int i, int br_i, int br_i1) {
	 	this.bp.addBreakPoint(source_index, i, br_i, br_i1);
		
	}
	public synchronized void addBreakPoint(int source_index, String i, int br_i, int br_i1) {
	 	this.bp.addBreakPoint(source_index, i, br_i, br_i1);
		
	}
	
	static Comparator<int[]> revcomp = new Comparator<int[]>() {

		@Override
		public int compare(int[] se1, int[] se2) {
			// TODO Auto-generated method stub
			if(se1[0]==se2[0]) return 0;
			return se1[0] < se2[0]? -1 : 1;
//			return -1*(Integer.compare(se1[0], se2[0]));
		}
		
		
	};
	
	static SortedMap<int[], String> copy(SortedMap<int[], SAMRecord> map) {
		TreeMap<int[], String> tm = new TreeMap<int[], String>(map.comparator());
		Iterator<Entry<int[], SAMRecord>> it = map.entrySet().iterator();
		while(it.hasNext()) {
			Entry<int[], SAMRecord> ent = it.next();
			SAMRecord sr = ent.getValue();
			int[] key = ent.getKey();
			tm.put(key,key[0]+","+key[1]+":"+
					sr.getReferenceName()+",neg:"+sr.getReadNegativeStrandFlag()+",supp:"+sr.isSecondaryOrSupplementary());
		}
		return tm;
	}
	static boolean PRINT_MAP=false;
	static boolean checkList(SortedMap<int[], SAMRecord> res, SAMRecord primary, 
			StringBuffers sb
			
			){
			String readname = primary.getReadName();
			String read_str = primary.getReadString();
			boolean primary_neg = primary.getReadNegativeStrandFlag();
		//	Iterator<int[]> it = res.keySet().iterator();
		//	int[] se_prev = it.next();
			int[] se_prev=null;
			boolean removed=false;
			String read_str1 = primary_neg ?  NeedlemanWunsch.reverseComplement(read_str) : read_str;
			int cnt=1;
			
			if(PRINT_MAP) {
				System.err.println(gson.toJson(res.keySet()));
				System.err.println(gson.toJson(copy(res).values()));
				if(false) {
					List<Boolean> strands=res.values().stream().map(t -> t.getReadNegativeStrandFlag()).collect(Collectors.toList());
					List<List<int[]>> strands1=res.values().stream().map(t -> 
					t.getAlignmentBlocks().stream().map(t1 -> stEnd(t1)).collect(Collectors.toList())).collect(Collectors.toList());
					List<int[]> strands2=res.values().stream().map(t -> stEnd1(t)).collect(Collectors.toList());
	
					System.err.println("h");
					}
			}
			byte[] bq = primary.getBaseQualities();
			Iterator<Entry<int[], SAMRecord>> it1 =  res.entrySet().iterator();
		
			//String chr = null;
			
			while(it1.hasNext()) {
				boolean first = se_prev==null;
				
				Entry<int[],SAMRecord> ent = it1.next();
				int[] se1 = ent.getKey();
				if(!first) {
					sb.append(";");
				
					int overlap = se_prev[1] - se1[0];
					sb.overlaps.append(overlap);
					SAMRecord sr1 = res.get(se_prev);
					SAMRecord sr2 = res.get(se1);
					String char1 = sr1.getReadNegativeStrandFlag() ? "-" : "+";
					String char2 = sr2.getReadNegativeStrandFlag() ? "-" : "+";
					boolean same =sr1.getReferenceName().equals(sr2.getReferenceName());
					String header = ">"+readname+",overlap:"+overlap+",cnt:"+cnt+","+se1[0]+","+se_prev[1]+","+char1+"," +char2+","+sr1.getReferenceName()+","+sr2.getReferenceName()+","+same+"\t";
					if(overlap>0) {
						String substr = read_str1.substring(se1[0]-1,se_prev[1]-1);
						if(Outputs.overlapOut!=null) {
							Outputs.overlapOut.print(header);
							Outputs.overlapOut.println(substr);
						}
						if(overlap>overlap_max) {
							removed=true;
						}
						//return null;

					}else if(overlap <0) {
						
						String substr = read_str1.substring(se_prev[1]-1, se1[0]-1);
						

						if(Outputs.joinOut!=null) {
							Outputs.joinOut.print(header);
							Outputs.joinOut.println(substr);
						}
						if(overlap < -1 *gap_max) {
							removed=true;
						}
					}else if(overlap<2 && overlap >-2) {
						String substr = read_str1.substring(se_prev[1]-6, se1[0]+4); //because its one based
						if(Outputs.noGap!=null) {
							Outputs.noGap.print(header);
							Outputs.noGap.println(substr);
						
						}
						
					}
					cnt++;
				}
				se_prev = se1;
				
				SAMRecord record = ent.getValue();
						sb.strand.append(record.getReadNegativeStrandFlag() ? '-': '+');
					
				sb.q_str.append(record.getMappingQuality());
				
				
			//	int[] se1 = ent.getKey();
				int qstart = se1[0];//record.getReadPositionAtReferencePosition(record.getAlignmentStart());
				int qend = se1[1];//record.getReadPositionAtReferencePosition(record.getAlignmentEnd());
				sb.read_pos.append(qstart+"_"+qend);
				
				//System.err.println("qstart-end "+qstart+" "+qend);
				sb.q_str1.append(TranscriptUtils1.median(bq, qstart, qend-qstart));
				
			 	sb.chrom.append(record.getReferenceName());
			 
			}
			
			
			
			
			return removed;
			
	}
	
	SortedMap<int[], SAMRecord> getList(List<SAMRecord> li){
		SAMRecord primary = li.get(0);
		if(primary.isSecondaryAlignment()) throw new RuntimeException("!");
		//sam1.clear();//seXM_004737643.2.clear();
		//sam1.addAll(sam_1);

		int read_length = primary.getReadLength();
		String read_str = primary.getReadString();
		boolean primary_neg = primary.getReadNegativeStrandFlag();
		TreeMap<int[], SAMRecord> res = new TreeMap<int[],SAMRecord>(revcomp);
	//	System.err.println(li.get(0).getReadName());
		//List<Boolean >res1 = new ArrayList<Boolean>();
		if(li.size()==1) {
			SAMRecord sr = li.get(0);
			res.put(new int[] {0, sr.getReadLength()}, sr);
			return res;
		}
		for(int k=0; k<li.size(); k++) {
			int[] stp = getStEnd(li.get(k), read_str,read_length, primary_neg);
			res.put(stp,(li.get(k)));
			
		}
		
		return res;
		
	}
	static int[] stEnd(AlignmentBlock ab) {
		int st = ab.getReadStart();
		int end = st+ab.getLength();
		return new int[] {st,end};
	}
	static int[] stEnd1(SAMRecord sr) {
		int st = sr.getAlignmentBlocks().get(0).getReadStart();
		AlignmentBlock ab = 
				sr.getAlignmentBlocks().get(sr.getAlignmentBlocks().size()-1);
		int end = ab.getReadStart()+ab.getLength();
		return new int[] {st,end};
	}
	public void empty(){
		int len =idents.size();
		//if(len>1) {
			//throw new RuntimeException(len+" !!");
	//	}
		for(int i=0; i<len; i++){
			 IdentityProfile1 idp=idents.get(i);
			// if(ViralTranscriptAnalysisCmd2.allowSuppAlignments) idp.commit();
			idp.clear();
		}
	}
	SortedSet<String> geneNames = new TreeSet<String>();
	
		
	static class StringBuffers{
		StringBuffer chrom = new StringBuffer();
		StringBuffer strand = new StringBuffer();//sam1.get(0).getReadNegativeStrandFlag() ? '-': '+';
		StringBuffer q_str = new StringBuffer();
		StringBuffer q_str1 = new StringBuffer();
		StringBuffer read_pos = new StringBuffer();
		StringBuffer overlaps = new StringBuffer();
		public void append(String string) {
			read_pos.append(";");
			chrom.append(";");//strand.append(";");
			q_str.append(";");q_str1.append(";");
			overlaps.append(";");
			// TODO Auto-generated method stub
			
		}	
	}
	

	/*SortedMap<Integer, Integer> currentStart = new TreeMap<Integer, Integer>();
	private synchronized void removeStart(int startPos){
		int  v1 = currentStart.get(startPos);
		if(v1==0) currentStart.remove(startPos);
		else currentStart.put(startPos, v1-1);
	}
	private synchronized void addStart(int startPos){
		Integer v = currentStart.get(startPos);
		currentStart.put(startPos, v==null ? 1 : v+1);
	}*/
	final Integer[] basesStart = new Integer[4];final Integer[] basesEnd = new Integer[4];
	final Integer[] readCounts = new Integer[] {0,0,0};
	static Alphabet alph = Alphabet.DNA();
	public void identity1(List<SAMRecord>sam_1,
			 boolean cluster_reads,     List<Sequence> genome, final int primaryIndex
			) {
		if(primaryIndex<0) {
			System.err.println("no primary index");
			return;
		}
		if(primaryIndex>0){
	System.err.println(" problem with primary index");			
return ;		
	//throw new RuntimeException("expected at zero");
		}
		//List<SAMRecord>sam_1 = new ArrayList<SAMRecord>();
		//sam_1.addAll(Arrays.asList(sam));
		if(sam_1.size()==0){
			return ;
		}
		// TODO Auto-generated method stub
		
		//final char strand2 = strand;
	//	final int start = sam.get(0).getStart();
	//	addStart(start);
	//	if(ViralTranscriptAnalysisCmd2.sorted && (!ViralTranscriptAnalysisCmd2.sequential || num_sources==1)  && this.all_clusters.l.size()> clearUpThreshold){ 
	//	 this.all_clusters.clearUpTo(currentStart.firstKey() -cleanUpMoreThan , o);
	//	}
		//Boolean backwardStrand =null;
		
	
		
		Runnable run = new Runnable(){
			public void run(){
				if(sam_1.size()==0) throw new RuntimeException("!!");
				boolean include = process(sam_1.get(primaryIndex));
				if(!include){
				
					return ;
				}
				
				
				String readSeq = sam_1.get(0).getReadString();//new Sequence(alph, sam_1.get(0).getReadString(), sam_1.get(0).getReadName());
///DONT THINK NEED TO SORT
			//	Collections.sort(sam_1 , comp_q);
				if(sam_1.get(0).isSecondaryOrSupplementary()) {
					throw new RuntimeException("first should be primary");
				}
		//		int source_index = (Integer) sam_1.get(0).getAttribute(SequenceUtils.src_tag);
				String rn= sam_1.get(0).getReadName();
				final IdentityProfile1 profile = get();
				int sze = sam_1.size();
				//List<SAMRecord>sam1 = new ArrayList<SAMRecord>();
		//	List<int[]>se = new ArrayList<int[]>();
				int limit = 1;
			//	inner: for(int j=0; sam_1.size()>0 && j<limit; j++){
					SAMRecord primary = sam_1.get(0);
					if(primary.isSecondaryAlignment()) throw new RuntimeException("!");
					//sam1.clear();//se.clear();
					//sam1.addAll(sam_1);
			
					int read_length = primary.getReadLength();
					String read_str = primary.getReadString();
			//	String read_strand = (String) primary.getAttribute(PolyAT.read_strand_tag);
					
					
					byte[] bq = primary.getBaseQualities();
				//	System.err.println(sam1.size()+" "+sze);
					
				//	boolean fusion = false;
					//boolean primary_neg = primary.getReadNegativeStrandFlag();
					//boolean multi = sam_1.size()>1;
					SortedMap<int[], SAMRecord> map = getList(sam_1);
					List<Integer> overlaps = new ArrayList<Integer>();
					StringBuffers sb = new StringBuffers();
					boolean remove = checkList(map, primary, sb);
					
//					System.err.println(gson.toJson(copy(map)));
/*					if(multi) {
						Comparator samComparator = new Comparator<SAMRecord>() {

							@Override
							public int compare(SAMRecord o1, SAMRecord o2) {
								int[] se1 = getStEnd(o1, read_str,read_length, primary_neg);
								int[] se2 = getStEnd(o2,read_str, read_length, primary_neg);
								return -1*(Integer.compare(se1[1], se2[1]));
							}
							
							
						};
					Collections.sort(sam_1, samComparator); // orders along the read
					}*/
					if(!remove) {
					
				
				
					int source_index=0;
					 profile.setName(rn, sb,source_index, map.size()>1);
				//	if(profile.coRefPositions.breaks.size()>0 && ! supp) throw new RuntimeException("this is not clear");
				//	if(supp && profile.suppl!=null && profile.suppl.size()>0) throw new RuntimeException("supps not empty");
					
					profile.identity1(genome,  readSeq, sam_1, source_index, cluster_reads);
					}
				//	profile.commit();
					
					
				//	removeStart(start);
				replace(profile);
			}

		

			private Object getStart(SAMRecord t) {
				// TODO Auto-generated method stub
				return null;
			}

			
		};
		if(executor==null) run.run();
		else executor.execute(run);
	

	}	
	
	public static void shutDownExecutor() {
		if(executor!=null) {
			waitOnThreads(10);
			executor.shutdown();
		}
		
	}
	


	private static boolean process(SAMRecord sam) {
		boolean reverseForward = ViralTranscriptAnalysisCmd2.reverseForward;
	//	int source_index = (Integer) sam.getAttribute(SequenceUtils.src_tag);
		//boolean FLAMES=ViralTranscriptAnalysisCmd2.barcodes==null  ? false: ViralTranscriptAnalysisCmd2.barcodes[source_index].FLAMES_BARCODES;

		boolean RNA = ViralTranscriptAnalysisCmd2.RNA;//[source_index];
		if(sam.isSecondaryOrSupplementary()){
			throw new RuntimeException();
		}
		int[] start_end_for = new int[2];
		int[] start_end_rev = new int[2];
		int[] start_for_pA = new int[2];
		int[] start_rev_pA = new int[2];
		
		boolean align_reverse = sam.getReadNegativeStrandFlag();

			String sa = sam.getReadString();
			if(ViralTranscriptAnalysisCmd2.illumina){
				sa = sam.getReadName();
			
			}else{
				if(align_reverse){  
					// this converts read back to original orientation
					sa = SequenceUtil.reverseComplement(sa);
				}
			}
			//Boolean forward_read=true;
			
//			System.err.println("RNA "+RNA[source_index]);
				//if(forward_read!=null){
					sam.setAttribute(PolyAT.read_strand_tag, "+");// forward_read ? "+" : "-");
				//}
		/*	
			if(ViralTranscriptAnalysisCmd2.barcodes!=null){
				if(forward_read==null) return false;
				
				if(RNA) throw new RuntimeException("not currently supporting RNA barcodes (but probably could)");
				try{
					int for_min;
				if(forward_read){
					for_min = ViralTranscriptAnalysisCmd2. barcodes[source_index].assign( sam,sa, true, start_end_for);
				
				}else{
					for_min = ViralTranscriptAnalysisCmd2.barcodes[source_index].assign( sam,sa, false, start_end_rev);
				}
				if( ViralTranscriptAnalysisCmd2.exclude_reads_without_barcode && for_min> Barcodes.tolerance_barcode){
					return false;
				}
				
			
				//if(forward_read!=null && forward_read1!=null){
					if(for_min <= Barcodes.tolerance_barcode) { //forward_read.equals(forward_read1)){
						int startpA;
						int startB;
						if(forward_read){ //reverse barcode at 3'end
							 startB =start_end_for[1];
							 startpA = start_for_pA[1];
						}else{
							 startB =start_end_rev[0];
							 startpA = start_rev_pA[0];
							
						}
						int diff = startpA-startB ;
					
//						System.err.println(forward_read+" "+startB+" "+startpA+" "+diff);
						if(diff< ViralTranscriptAnalysisCmd2.max_umi && diff>ViralTranscriptAnalysisCmd2.min_umi && !FLAMES){
							String umi;
							if(forward_read){
								
								//String barcode = (String)sam.getAttribute(Barcodes.barcode_forward_tag);
								//System.err.println(barcode);
								 int read_len = sa.length();
								umi = SequenceUtil.reverseComplement(sa.substring(Math.max(0,read_len - startpA), Math.min(read_len,read_len - startB)));
								//String umi1 = SequenceUtil.reverseComplement(sa.substring(read_len - startpA, read_len - startB+1));
							//	System.err.println(umi);
								//System.err.println(umi1);
							//	System.err.println("");
							}else{
								umi =   sa.substring(startB, startpA);
							}
						//	System.err.println(forward_read+" "+startB+" "+startpA+" "+diff+" "+umi);

							sam.setAttribute(Barcodes.umi_tag, umi);
						}
						
						sam.setAttribute(Barcodes.barcode_confidence_tag, "high");
					}else{
						sam.setAttribute(Barcodes.barcode_confidence_tag, "low");
					}
			//	}else{
				//	sam.setAttribute(Barcodes.barcode_confidence_tag, "NA");
			//	}
				
				}catch(Exception exc){
					exc.printStackTrace();
				}
			}
			*/
		
			return true;
	}

	
 
}
