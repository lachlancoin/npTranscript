package npTranscript.cluster;

import java.io.IOException;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.SequenceUtil;
import japsa.bio.np.barcode.SWGAlignment;
import japsa.seq.Alphabet;
import japsa.seq.Sequence;


public class TranscriptUtils {
	
	public static int break_thresh = 1000;
	public static int endThresh = 100;
	public static int startThresh = 100;
	public static int extra_threshold = 500;
	public static int extra_threshold1 = 50;
	public static int extra_threshold2 = 10;

	public static boolean coronavirus = true;
	static int round(int pos, int round) {
		int res = (int) Math.floor((double) (pos) / (double) round);
		return res;
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
	
	
	
	/**
	 * Get the identity between a read sequence from a sam and a reference sequence
	 * 
	 * @param refSeq
	 * @param sam
	 * @return
	 */
	public static void identity1(Sequence refSeq, Sequence fivePrimeRefSeq, Sequence threePrimeRefSeq,   Sequence readSeq, SAMRecord sam, IdentityProfile1 profile, 
			int source_index, boolean cluster_reads, int seqlen) throws NumberFormatException{
		int readPos = 0;// start from 0
		//Outputs output = profile.o;
		int refPos = sam.getAlignmentStart() - 1;// convert to 0-based index
		String id = sam.getReadName();
		
		profile.newRead(source_index);
		for (final CigarElement e : sam.getCigar().getCigarElements()) {
			final int length = e.getLength();

			switch (e.getOperator()) {
			case H:
				// nothing todo
			//	profile.readClipped += length;
				break; // ignore hard clips
			case P:
			//	profile.readClipped += length;
				// pad is a kind of hard clipped ??
				break; // ignore pads
			case S:
				// advance on the reference
				//profile.readClipped += length;
				readPos += length;
				break; // soft clip read bases
			case N:
				// System.err.println(length);
				refPos += length;
				/*for (int i = 0; i < length && refPos + i < refSeq.length(); i++) {
					profile.refClipped[refPos + i] += 1;
				}*/
				// profile.refClipped += length;
				break; // reference skip

			case D:// deletion
				refPos += length;
				profile.refBase += length;
				for (int i = 0; i < length && refPos + i < refSeq.length(); i++) {
				//	profile.baseDel[refPos + i] += 1;
					profile.addRefPositions(refPos + i, false);
				}
				//profile.numDel++;
				break;

			case I:
				readPos += length;
				profile.readBase += length;

			//	profile.baseIns[refPos] += length;
			//	profile.numIns++;
				break;
			case M:
				for (int i = 0; i < length && refPos + i < refSeq.length(); i++) {
					
					if (refSeq.getBase(refPos + i) == readSeq.getBase(readPos + i)){
				//		profile.match[refPos + i]++;
						profile.addRefPositions(refPos + i, true);
					}
					else{
					//	profile.mismatch[refPos + i]++;
						profile.addRefPositions(refPos + i, false);
					}
				}
				profile.readBase += length;
				profile.refBase += length;

				readPos += length;
				refPos += length;
				break;

			case EQ:
				readPos += length;
				refPos += length;
				for (int i = 0; i < length && (refPos + i) < refSeq.length(); i++) {
					profile.addRefPositions(refPos+i, true);
				}
				profile.readBase += length;
				profile.refBase += length;
				//profile.match[refPos] += length;
				break;

			case X:
				readPos += length;
				refPos += length;

				profile.readBase += length;
				profile.refBase += length;
				for (int i = 0; i < length && (refPos + i) < refSeq.length(); i++) {
					profile.addRefPositions(refPos+i, false);
				}
				//profile.addRefPositions(refPos);
			//	profile.mismatch[refPos] += length;
				break;
			default:
				throw new IllegalStateException("Case statement didn't deal with cigar op: " + e.getOperator());
			}// case
		} // for
		try{
			int maxl = 100;
			int tol=5;
			int st_r = sam.getReadPositionAtReferencePosition(sam.getAlignmentStart());
			int end_r = sam.getReadPositionAtReferencePosition(sam.getAlignmentEnd());
			String baseQ = sam.getBaseQualityString();
			char strand = sam.getReadNegativeStrandFlag() ? '-': '+';
			SWGAlignment align_5prime = null;
			SWGAlignment align_3prime = null;
			SWGAlignment align_3primeRev = null;

			int offset_3prime =0;
			int polyAlen = TranscriptUtils.polyAlen(refSeq);
			if(TranscriptUtils.coronavirus && st_r > extra_threshold1 && st_r < 1000 && sam.getAlignmentStart()> 100 ){
				//good candidate for a missed 5' alignment
				 align_5prime = SWGAlignment.align(readSeq.subSequence(0, st_r), fivePrimeRefSeq);
			}
			int diff_r = readSeq.length() -end_r;
			int seqlen1 = refSeq.length()-polyAlen;
			if(TranscriptUtils.coronavirus && diff_r > extra_threshold2 && diff_r < 1000 && sam.getAlignmentEnd()< seqlen1- 100 ){
				//good candidate for a missed 3' alignment
				 align_3prime = SWGAlignment.align(readSeq.subSequence(end_r+1,readSeq.length()), threePrimeRefSeq);
				// align_3primeRev = SWGAlignment.align(TranscriptUtils.revCompl(readSeq.subSequence(end_r+1,readSeq.length())), threePrimeRefSeq);

				 offset_3prime = refSeq.length()-threePrimeRefSeq.length();
			}
			//String flanking = start_flank1+"\t"+start_flank2+"\t"+end_flank1+"\t"+end_flank2;
			//SAMRecord sam, String id, boolean cluster_reads, int  readLength, Sequence refSeq, int src_index , Sequence readSeq, String baseQ, 
			//int start_read, int end_read, char strand, SWGAlignment align5prime
			
			boolean splice = profile.processRefPositions(sam, id, cluster_reads, 
					refSeq, source_index, readSeq,baseQ, st_r, end_r, strand, align_5prime, align_3prime,align_3primeRev, offset_3prime, polyAlen);
			//String ID = profile.clusterID
			if(!splice && Outputs.writeUnSplicedFastq){
				
				String start_flank1 = refSeq.subSequence(Math.max(0, sam.getAlignmentStart()-10),sam.getAlignmentStart()).toString();
				String start_flank2 = refSeq.subSequence(sam.getAlignmentStart(),Math.min(sam.getAlignmentStart()+10, refSeq.length())).toString();
				String end_flank1 = refSeq.subSequence(Math.max(0, sam.getAlignmentEnd()-10),sam.getAlignmentEnd()).toString();
				String end_flank2 = refSeq.subSequence(sam.getAlignmentEnd(),Math.min(sam.getAlignmentEnd()+10,refSeq.length())).toString();
				
				int len = sam.getAlignmentEnd() - sam.getAlignmentStart();
				readSeq.setDesc(profile.clusterID[0]+" "+profile.clusterID[1]+" "+st_r+" "+end_r+" "+readSeq.length()+" "+sam.getAlignmentStart()+" "+sam.getAlignmentEnd()+" "+len+" "+start_flank1+" "+start_flank2+" "+end_flank1+" "+end_flank2);
				profile.o.writeUnspliced(readSeq, baseQ, sam.getReadNegativeStrandFlag(), source_index);
			}
		//	System.err.println(refSeq.length());
			if(!splice && (sam.getAlignmentStart()>100 || sam.getAlignmentEnd()<(refSeq.length()-100))){
			String desc = ";start="+sam.getAlignmentStart()+";end="+sam.getAlignmentEnd()+";full_length="+readSeq.length()+";strand="+strand;
			
			if(st_r >extra_threshold  ){
				Sequence leftseq = readSeq.subSequence(0, st_r);
				String baseQL = baseQ.equals("*") ?  baseQ : baseQ.substring(0, st_r);

				if(sam.getReadNegativeStrandFlag()) {
					leftseq =TranscriptUtils.revCompl(leftseq); 
							
					baseQL = (new StringBuilder(baseQ)).reverse().toString();
				}
			
				leftseq.setName(readSeq.getName()+".L");
				
				leftseq.setDesc("len="+leftseq.length()+desc);//+";mtch5="+mtch_5+";mtch_3="+mtch_3);
				profile.o.writeLeft(leftseq,baseQL,sam.getReadNegativeStrandFlag(), source_index);// (double) Math.max(mtch_3, mtch_5)> 0.7 * (double)leftseq1.length());
				
			}
			if(readSeq.length() -  end_r >extra_threshold){
				Sequence rightseq = readSeq.subSequence(end_r, readSeq.length());
				Sequence spanning1 = refSeq.subSequence(Math.max(0, sam.getAlignmentEnd()-10),sam.getAlignmentEnd());
				Sequence spanning2 = refSeq.subSequence(sam.getAlignmentEnd(),Math.min(refSeq.length(), sam.getAlignmentEnd()+10));
				
				 rightseq.setName(readSeq.getName()+".R");
					String baseQR = baseQ.equals("*") ?  baseQ : baseQ.substring(end_r, readSeq.length());
				 
				if(sam.getReadNegativeStrandFlag()) {
					rightseq = TranscriptUtils.revCompl(rightseq);
					baseQR = (new StringBuilder(baseQ)).reverse().toString();
				
				}
				 align_5prime = SWGAlignment.align(rightseq, refSeq);
				 
				
				//int mtch_5  = TranscriptUtils.checkAlign ?  checkAlignmentIsNovel(rightseq1, refSeq5prime, "right 5") : 0;
				//int mtch_3 = TranscriptUtils.checkAlign ? checkAlignmentIsNovel(rightseq1, refSeq3prime, "right 3"): 0;
				rightseq.setDesc("len="+rightseq.length()+desc+" "+spanning1.toString()+" "+spanning2.toString());//+";mtch5="+mtch_5+";mtch_3="+mtch_3);
				profile.o.writeLeft(rightseq,baseQR, sam.getReadNegativeStrandFlag(), source_index);///,  (double) Math.max(mtch_3, mtch_5)> 0.7 * (double)rightseq1.length());
			}
			
		}
		}catch(IOException exc){
			exc.printStackTrace();
		}

	}
	
	private static Sequence revCompl(Sequence leftseq) {
		String sequence = SequenceUtil.reverseComplement(leftseq.toString());
		return new Sequence(Alphabet.DNA(), sequence.toCharArray(), leftseq.getName());
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

	private static int checkAlignmentIsNovel(Sequence leftseq, Sequence refSeq, String type) {
	
			SWGAlignment align = SWGAlignment.align(leftseq, refSeq);
			return align.getIdentity();
		
			
		
	}
	
	

}