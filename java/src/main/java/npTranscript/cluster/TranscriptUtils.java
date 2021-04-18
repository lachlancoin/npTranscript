package npTranscript.cluster;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.SequenceUtil;
import japsa.bio.np.barcode.SWGAlignment;
import japsa.seq.Alphabet;
import japsa.seq.Sequence;


public class TranscriptUtils {
	
	public static boolean coronavirus = true;
	static int round(int pos, int round) {
		int res = (int) Math.floor((double) (pos) / (double) round);
		return res;
	}
	
	static int round(int pos, int round, boolean left) {
		int res = (int) Math.floor((double) (pos) / (double) round);
		return left ? res : res+1;
	}
	
//	static int[] printedLines =new int[] {0,0,0,0};
	public static int overlap(int st1, int end1, int st2, int end2){
		int minlen = Math.min(end1-st1, end2 - st2);
		int overlap = Math.min(end1-st2, end2 - st1);
		return Math.max(0, Math.min(minlen, overlap) +1);
	}
	public static int union(int st1, int end1, int st2, int end2, int overlap){
		int len1 = end1 - st1;
		int len2 = end2 - st2;
		if(overlap<=0) {
			return len1+len2;
		}else{
			int union =  Math.max(end1-st2, end2 - st1);
			int maxlen = Math.max(len1, len2);
			return Math.max(union, maxlen)+1 - overlap;
		}
	}
	
	static String getString(int[] c) {
		StringBuffer sb = new StringBuffer(c[0]+"");
		for(int i=1; i<c.length;i++) {
			sb.append("\t");sb.append(c[i]);
		}
		return sb.toString();
	}
	static String getString(double[] c) {
		String formatstr = "%5.3g";
		StringBuffer sb = new StringBuffer(String.format(formatstr,c[0]).trim()+"");
		for(int i=1; i<c.length;i++) {
			sb.append("\t");sb.append(String.format(formatstr, c[i]).trim());
		}
		return sb.toString();
	}
	 static String getString(String string, int num_sources2, boolean print_index) {
		 StringBuffer sb = new StringBuffer(string);
		 if(print_index) sb.append(0);
			for(int i=1; i<num_sources2;i++) {
				sb.append("\t");sb.append(string);
				if(print_index)sb.append(i);
			}
			String res = sb.toString();
			return res;
	}
//	static double round2 = 100;
	
	static String[] nmes = "5_3:5_no3:no5_3:no5_no3".split(":");
	public static boolean checkAlign = true;
	public  static boolean checkNegStrand = false;
	
	
	
	
	
	//public static boolean findPolyA = false;
	
	
	public static int endThresh = 100;
	public static int startThresh = 100;

	public static Sequence rev(Sequence leftseq) {
		String sequence = (leftseq.toString());
		return new Sequence(Alphabet.DNA(), (new StringBuilder(sequence)).reverse().toString().toCharArray(), leftseq.getName());
	}
	public static Sequence compl(Sequence leftseq) {
		String sequence = SequenceUtil.reverseComplement(leftseq.toString());
		return new Sequence(Alphabet.DNA(), (new StringBuilder(sequence)).reverse().toString().toCharArray(), leftseq.getName());
	}

	public static Sequence revCompl(Sequence leftseq) {
		String sequence = SequenceUtil.reverseComplement(leftseq.toString());
	   Sequence res = new Sequence(Alphabet.DNA(), sequence.toCharArray(), leftseq.getName());
	  res.setDesc(leftseq.getDesc());
	   return res;
	}
	
	public static void reverseArray(byte[] intArray) {
		 int size = intArray.length;
	        int i, k;
	        byte temp; 
	        for (i = 0; i < size / 2; i++) { 
	            temp = intArray[i]; 
	            intArray[i] = intArray[size - i - 1]; 
	            intArray[size - i - 1] = temp; 
	        } 
	 } 
	
	public static int polyAlen(Sequence refSeq){
		int seqlen = refSeq.length();
		char[] last10bp = refSeq.subSequence(seqlen-10, seqlen).charSequence();
		boolean polyA = true;
		for(int i=0; i<last10bp.length; i++){
			if(last10bp[i]!='A') polyA = false;
			
		}
		if(polyA){
			int i =1;
			while(refSeq.charAt(seqlen-i)=='A'){
				i++;
			}
			return i-1;
		}
		return 0;
	}
	
			;
	public static void count(String read, boolean start, Integer[] bases, int len, String chars){
		
		//String read = revCompl ? TranscriptUtils.
		int seqlen = read.length();
		//String chars_ = revCompl? chars_revC : chars;
//		boolean start = revCompl ? !start_ : start_;
		
		for(int i=0; i<len; i++){ // should it start at 0
			int i1  = start ? i : seqlen-i-1;
			char c = read.charAt(i1);
			
			bases[chars.indexOf(c) ]++;
		}
		
	}

	private static String getString(SWGAlignment align_5prime, boolean neg) {
		 int st1 = align_5prime.getStart1();
		 int st2 = align_5prime.getStart2();
		 int len1 = align_5prime.getSequence1().length - align_5prime.getGaps1();
		 int len2 = align_5prime.getSequence2().length -align_5prime.getGaps2();
		 int end1 = st1 + len1;
		 int end2 = st2 + len2;
	//	 System.err.println(align)
		 char strand = neg ? '-' : '+';
		 return strand+","+align_5prime.getScore()+","+align_5prime.getIdentity()+" "+st1+","+end1+" "+st2+","+end2;
	}

	

	public static boolean writeAnnotP = true;



	//public static String bedChr= null;
	
	 static  void getStartEnd(SWGAlignment align, Integer[] seq1, Integer[] seq2,  int offset1, int offset2,  boolean negStrand1) {
		int len1 = align.getOriginalSequence1().length();
		seq1[0] = align.getStart1();
		seq1[1] = align.getStart1()+align.getSequence1().length - align.getGaps1();
		seq2[1] = align.getIdentity();
		seq2[0] = align.getLength();
		seq1[2] = align.getStart2();
		seq1[3] = align.getStart2()+align.getSequence2().length - align.getGaps2();
		//seq2[2] = align.getIdentity();
		//seq2[3] = align.getLength();
		if(negStrand1){
			int a = seq1[1];
			seq1[1] = len1 - seq1[0];
			seq1[0] = len1 - a;
		}
		seq1[0]= seq1[0]+offset1;
		seq1[1] = seq1[1]+offset1;
		seq1[2]= seq1[2]+offset2;
		seq1[3] = seq1[3]+offset2;
	}

	 /* flips read and quality without chaning the flag */
	public static void flip(SAMRecord sam, boolean switchFlag) {
		if(true) throw new RuntimeException ("!!");
		String sa = sam.getReadString();
		byte[]phredQs = sam.getBaseQualities();
		byte[] bases = sam.getReadBases();
		TranscriptUtils.reverseArray(phredQs);
		TranscriptUtils.reverseArray(bases);
		sam.setReadBases(bases);
		sam.setBaseQualities(phredQs);
		sa = SequenceUtil.reverseComplement(sa);
		sam.setReadString(sa);
		if(switchFlag){
			sam.setReadNegativeStrandFlag(!sam.getReadNegativeStrandFlag());
		}
	}
	
	

}