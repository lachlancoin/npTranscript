package npTranscript.cluster;

import java.io.IOException;
import java.util.Arrays;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.SequenceUtil;
import japsa.bio.np.barcode.SWGAlignment;
import japsa.seq.Alphabet;
import japsa.seq.Sequence;


public class TranscriptUtils {
	
	public static double qual_thresh = 20.0D;

	
	public static int break_thresh = 1000;
	public static int endThresh = 100;
	public static int startThresh = 100;
	public static int extra_threshold = 500;
	public static int extra_threshold1 = 50;
	public static int extra_threshold2 = 20;
	public static boolean attempt5rescue = false;
	public static boolean attempt3rescue = false;

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
	
	
	
	
	public static Sequence polyA = new Sequence(Alphabet.DNA(), "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA".toCharArray(), "polyA");
	public  static boolean tryComplementOnExtra = false;
	public static boolean reAlignExtra = false;
	//public static boolean findPolyA = false;
	/**
	 * Get the identity between a read sequence from a sam and a reference sequence
	 * 
	 * @param refSeq
	 * @param sam
	 * @return
	 */
	public static void identity1(Sequence refSeq, Sequence fivePrimeRefSeq, Sequence threePrimeRefSeq,   Sequence readSeq, SAMRecord sam, IdentityProfile1 profile, 
			int source_index, boolean cluster_reads, int seqlen, String poolID) throws NumberFormatException{
		int readPos = 0;// start from 0
		//Outputs output = profile.o;
		int refPos = sam.getAlignmentStart() - 1;// convert to 0-based index
		String id = sam.getReadName();
		String chrom = refSeq.getName();
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
			int diff_r = readSeq.length() -end_r;
			
			
			
			String baseQ = sam.getBaseQualityString();
			byte[]phredQs = sam.getBaseQualities();
			char strand = sam.getReadNegativeStrandFlag() ? '-': '+';
			boolean forward = !sam.getReadNegativeStrandFlag();
			SWGAlignment align_5prime = null;
			SWGAlignment align_3prime = null;
			SWGAlignment align_3primeRev = null;

			int offset_3prime =0;
			int polyAlen = TranscriptUtils.polyAlen(refSeq);
			if(TranscriptUtils.coronavirus && st_r > extra_threshold1 && attempt5rescue && st_r < 1000 && sam.getAlignmentStart()> 100 ){
				//good candidate for a missed 5' alignment
				try{
				 align_5prime = SWGAlignment.align(readSeq.subSequence(0, st_r), fivePrimeRefSeq);
				}catch(Exception exc){
					System.err.println("warning could not attemp 5'rescue "+ readSeq.subSequence(0, st_r));
				}
			}
			
			int seqlen1 = refSeq.length()-polyAlen;
			if(TranscriptUtils.coronavirus && diff_r > extra_threshold2 && attempt3rescue &&  diff_r < 1000 && sam.getAlignmentEnd()< seqlen1- 100 ){
				//good candidate for a missed 3' alignment
				try{
				 align_3prime = SWGAlignment.align(readSeq.subSequence(end_r+1,readSeq.length()), threePrimeRefSeq);
				// align_3primeRev = SWGAlignment.align(TranscriptUtils.revCompl(readSeq.subSequence(end_r+1,readSeq.length())), threePrimeRefSeq);

				 offset_3prime = refSeq.length()-threePrimeRefSeq.length();
				}catch(Exception exc){
					System.err.println("warning could not attempt 3'rescue "+ readSeq.subSequence(end_r+1,readSeq.length()));
				}
			}
			Sequence leftseq = null;  Sequence rightseq = null;
		 double phredQL = st_r < 20 ? 0 : npTranscript.run.ViralChimericReadsAnalysisCmd.median(phredQs, 0,st_r);
		// double phredQ = npTranscript.run.ViralChimericReadsAnalysisCmd.median(phredQs, st_r, end_r - st_r);
		 double phredQR = diff_r <20 ? 0 : npTranscript.run.ViralChimericReadsAnalysisCmd.median(phredQs,  end_r , diff_r);
		 //System.err.println(phredQL+" "+phredQ+" "+phredQR);
			if(st_r> extra_threshold && Outputs.writePolyA && !Double.isNaN(phredQL) && phredQL >= qual_thresh){
				leftseq = readSeq.subSequence(0,st_r);
				try{
				SWGAlignment polyAlign =  SWGAlignment.align(leftseq, polyA);
				if(polyAlign.getIdentity() > 0.9 * polyAlign.getLength()  && polyAlign.getLength()>15){
					int st = polyAlign.getStart1();
					int end = st + polyAlign.getLength() - polyAlign.getGaps1();
					if(st> 10){
						String nme = readSeq.getName();
					//	readSeq.setName(nme+".L");
						profile.o.writePolyA(readSeq.subSequence(0, end),  nme+".L",
								baseQ.length()==1 ? baseQ: baseQ.substring(0,end), sam.getReadNegativeStrandFlag(),source_index);
					//	readSeq.setName(nme+".R");
						profile.o.writePolyA(readSeq.subSequence( end+1, readSeq.length()),  nme+".R",
								baseQ.length()==1 ? baseQ: baseQ.substring(end+1, readSeq.length()), sam.getReadNegativeStrandFlag(),source_index);
						System.err.println("internal polyA 5´ ");
						return;
					}
				}
				}catch(Exception exc){
					System.err.println("warning could not do polyA alignment");
				}

			}
			if(readSeq.length() -  end_r >extra_threshold && Outputs.writePolyA &&  !Double.isNaN(phredQR) && phredQR >= qual_thresh){
				rightseq = readSeq.subSequence(end_r+1, readSeq.length());
				try{
				SWGAlignment polyAlign =  SWGAlignment.align(rightseq, polyA);
				if(polyAlign.getIdentity() > 0.9 * polyAlign.getLength()  && polyAlign.getLength()>10 ){
					 int st = polyAlign.getStart1() + end_r;
					 int end = st + polyAlign.getLength() - polyAlign.getGaps1();
					 if(end< readSeq.length()-10){
						 String nme = readSeq.getName();
						
							profile.o.writePolyA(readSeq.subSequence(0, end), nme+".L",
									baseQ.length()==1 ? baseQ: baseQ.substring(0,end), sam.getReadNegativeStrandFlag(),source_index);
							profile.o.writePolyA(readSeq.subSequence( end, readSeq.length()),  nme+".R",
									baseQ.length()==1 ? baseQ: baseQ.substring(end, readSeq.length()), sam.getReadNegativeStrandFlag(),source_index);
							System.err.println("internal polyA 3'");
					return;
					 }
				}
				}catch(Exception exc){
					System.err.println("warning could not do polyA alignment");
				}
			}
			 
			String secondKey= profile.processRefPositions(sam, id, cluster_reads, 
					refSeq, source_index, readSeq,baseQ, st_r, end_r, strand, align_5prime, align_3prime,align_3primeRev, offset_3prime, polyAlen, poolID);
			//String ID = profile.clusterID
			diff_r = readSeq.length() - profile.readEn;
			seq11[0]=  profile.readSt; seq11[1] = profile.readEn; 
			seq11[2] = profile.startPos; seq11[3] = profile.endPos;
			
		//	System.err.println(refSeq.length());
			//check read not mapping to start or end of reference (?)
			if( (st_r>extra_threshold || end_r>extra_threshold)){
			//String desc = ";start="+sam.getAlignmentStart()+";end="+sam.getAlignmentEnd()+";full_length="+readSeq.length()+";strand="+strand;
				
			if(st_r >extra_threshold  && st_r > diff_r  && !Double.isNaN(phredQL) && phredQL >= qual_thresh){
				if(leftseq==null) leftseq = readSeq.subSequence(0,st_r);
				//if(leftseq==null) leftseq = readSeq.subSequence(0, st_r);
				String baseQL = baseQ.equals("*") ?  baseQ : baseQ.substring(0, st_r);

				if(sam.getReadNegativeStrandFlag()) {
					leftseq =TranscriptUtils.revCompl(leftseq); 
					leftseq.setDesc(leftseq.getDesc()+" rev "+chrom);		
					baseQL = (new StringBuilder(baseQ)).reverse().toString();
				}
				else{
					leftseq.setDesc(leftseq.getDesc()+" fwd "+chrom);		
				}
			
				leftseq.setName(readSeq.getName());
				StringBuffer desc = new StringBuffer();
				desc.append("L "+secondKey+" "+st_r+" "+getString(seq11));
			
				if(reAlignExtra){
					Sequence refSeq1 = refSeq.length() < 30000  ? refSeq : 
							refSeq.subSequence(profile.startPos - 10000, profile.endPos+1000);
				align_5prime = SWGAlignment.align(leftseq, refSeq1);
				 TranscriptUtils.getStartEnd(align_5prime, seq1, seq2, 0, 0, sam.getReadNegativeStrandFlag());
				String secondKey1 =  profile.all_clusters.annot.nextUpstream(seq1[2], profile.chrom_index, forward)+";"+profile.all_clusters.annot.nextDownstream(seq1[3], profile.chrom_index,forward);
				desc.append(" "+String.format("%5.3g",(double)seq2[1]/(double) seq2[0]).trim()+" "+secondKey1+" "+getString(seq1)+";"+getString(seq2));
				if(tryComplementOnExtra){
					align_5prime = SWGAlignment.align(TranscriptUtils.revCompl(leftseq), refSeq);
					 TranscriptUtils.getStartEnd(align_5prime, seq1, seq2, 0, 0, !sam.getReadNegativeStrandFlag());
						 secondKey1 =  profile.all_clusters.annot.nextUpstream(seq1[2], profile.chrom_index,forward)+";"+profile.all_clusters.annot.nextDownstream(seq1[3], profile.chrom_index,forward);

					 desc.append(" "+secondKey1+" "+String.format("%5.3g",(double)seq2[1]/(double) seq2[0]).trim()+" "+getString(seq1)+";"+getString(seq2));
						
				}
				}
				
				leftseq.setDesc(desc.toString());

				//leftseq.setDesc("len="+leftseq.length()+desc);//+";mtch5="+mtch_5+";mtch_3="+mtch_3);
				profile.o.writeLeft(leftseq,baseQL,sam.getReadNegativeStrandFlag(), source_index);// (double) Math.max(mtch_3, mtch_5)> 0.7 * (double)leftseq1.length());
				
			}
			if(diff_r >extra_threshold && diff_r> st_r && !Double.isNaN(phredQR) && phredQR >= qual_thresh){
				if(rightseq==null)			rightseq = readSeq.subSequence(end_r+1, readSeq.length());
			//	Sequence rightseq = readSeq.subSequence(end_r, readSeq.length());
			///	Sequence spanning1 = refSeq.subSequence(Math.max(0, sam.getAlignmentEnd()-10),sam.getAlignmentEnd());
			//	Sequence spanning2 = refSeq.subSequence(sam.getAlignmentEnd(),Math.min(refSeq.length(), sam.getAlignmentEnd()+10));
				
				 rightseq.setName(readSeq.getName());
					String baseQR = baseQ.equals("*") ?  baseQ : baseQ.substring(end_r, readSeq.length());
				 
				if(sam.getReadNegativeStrandFlag()) {
					rightseq = TranscriptUtils.revCompl(rightseq);
					baseQR = (new StringBuilder(baseQ)).reverse().toString();
				
				}
				StringBuffer desc = new StringBuffer();
				desc.append("R "+secondKey+" "+diff_r+" "+getString(seq11));
				if(reAlignExtra){
					
					Sequence refSeq1 = refSeq.length() < 30000  ? refSeq : 
							refSeq.subSequence(profile.startPos - 10000, profile.endPos+1000);
					try{
				 align_3prime = SWGAlignment.align(rightseq, refSeq1);
					
				 TranscriptUtils.getStartEnd(align_3prime, seq1, seq2, end_r, 0, sam.getReadNegativeStrandFlag());
					String secondKey1 =  profile.all_clusters.annot.nextUpstream(seq1[2], profile.chrom_index,forward)+";"+profile.all_clusters.annot.nextDownstream(seq1[3], profile.chrom_index,forward);

				 desc.append(" "+String.format("%5.3g",(double)seq2[1]/(double) seq2[0]).trim()+" "+secondKey1+" "+getString(seq1)+";"+getString(seq2));
				 if(tryComplementOnExtra){
						align_3prime = SWGAlignment.align(TranscriptUtils.revCompl(rightseq), refSeq);
						 TranscriptUtils.getStartEnd(align_3prime, seq1, seq2, end_r, 0, !sam.getReadNegativeStrandFlag());
							 secondKey1 =  profile.all_clusters.annot.nextUpstream(seq1[2], profile.chrom_index,forward)+";"+profile.all_clusters.annot.nextDownstream(seq1[3], profile.chrom_index,forward);

						 desc.append(" "+secondKey1+" "+String.format("%5.3g",(double)seq2[1]/(double) seq2[0]).trim()+" "+getString(seq1)+";"+getString(seq2));
							
					}
					}catch(Exception exc){
						System.err.println("warning 3'alignment unsuccessful");
					}
				}
				 rightseq.setDesc(desc.toString());
				 profile.o.writeLeft(rightseq,baseQR, sam.getReadNegativeStrandFlag(), source_index);///,  (double) Math.max(mtch_3, mtch_5)> 0.7 * (double)rightseq1.length());
			}
			
		}
		}catch(IOException exc){
			exc.printStackTrace();
		}

	}
	
	private static String getString(Integer[] seq12) {
		return CigarHash2.getString(Arrays.asList(seq12));
	}

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

	static Integer[] seq1 = new Integer[4];
	static Integer[] seq2 = new  Integer[2];
	static Integer[] seq11 = new Integer[4];
	static Integer[] seq21 = new  Integer[2];


	public static boolean writeAnnotP = true;
	
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
	
	

}