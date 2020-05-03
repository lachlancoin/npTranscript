package npTranscript.cluster;

import java.io.IOException;
import java.util.HashSet;
import java.util.Set;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.SAMRecord;
import japsa.seq.Sequence;


public class TranscriptUtils {
	
	public static int break_thresh = 10;
	public static int endThresh = 100;
	public static int startThresh = 100;
	public static int extra_threshold = 500;
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
	
	
	
	/**
	 * Get the identity between a read sequence from a sam and a reference sequence
	 * 
	 * @param refSeq
	 * @param sam
	 * @return
	 */
	public static void identity1(Sequence refSeq, Sequence readSeq, SAMRecord sam, IdentityProfile1 profile, 
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
			
			//sam.
			//System.err.println(st+" "+end+" "+len);
		//	System.err.println(sam.getUnclippedStart()+" "+sam.getUnclippedEnd()+" "+len);
		int st_r = sam.getReadPositionAtReferencePosition(sam.getAlignmentStart());
		int end_r = sam.getReadPositionAtReferencePosition(sam.getAlignmentEnd());
		//System.err.println(st_r + " "+end_r);
		//readSeq.setDesc("len="+readSeq.length());
		if(st_r >extra_threshold || (readSeq.length() -  end_r )>extra_threshold ){
			String desc = ";start="+sam.getAlignmentStart()+";end="+sam.getAlignmentEnd()+";full_length="+readSeq.length();

			Sequence leftseq = readSeq.subSequence(0, st_r);
			leftseq.setName(readSeq.getName()+".L");
			leftseq.setDesc("len="+leftseq.length()+desc);
			if(leftseq.length() > extra_threshold) profile.o.writeLeft(leftseq,source_index);

			Sequence rightseq = readSeq.subSequence(end_r, readSeq.length());
			rightseq.setName(readSeq.getName()+".R");
			rightseq.setDesc("len="+rightseq.length()+desc);
			if(rightseq.length() > extra_threshold) profile.o.writeLeft(rightseq,source_index);
			
			/*Sequence midseq = readSeq.subSequence(st_r, end_r);
			midseq.setName(readSeq.getName()+".M");
			midseq.setDesc("len="+midseq.length()+desc);
			profile.o.writeLeft(midseq,source_index);
			*/
			//profile.o.writeLeft(readSeq);
			
			//subseq.setDesc("st="+st_r);
			//System.err.println(subseq.toString());
		}
		
		profile.processRefPositions(sam.getAlignmentStart(), sam.getAlignmentEnd(), id, cluster_reads, readSeq.length(), refSeq.length(), source_index, readSeq,st_r, end_r);
		}catch(IOException exc){
			exc.printStackTrace();
		}

	}
	
	

}