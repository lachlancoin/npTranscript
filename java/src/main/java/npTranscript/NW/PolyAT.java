package npTranscript.NW;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.zip.GZIPInputStream;

import edlib.EdlibAlignConfig;
import edlib.EdlibAlignResult;
import edlib.EdlibAlignResult.ByValue;
import edlib.EdlibLibrary;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.SequenceUtil;
import japsa.seq.Sequence;
import npTranscript.cluster.TranscriptUtils;

public class PolyAT {
	private static boolean verbose=false;
	private static  int bc_len_AT = 6; //perhaps should be 15. seems shorter for RNA
	public static int edit_thresh_AT=1; // how many mismatches
	public static int max_dist_to_end_AT=20;  //how far from ends to search  150 for cDNA 100 for dRNA
	
	public static int dist_to_ignore_AT = 0;
	
	public static String polyAT_forward_tag = "AA";
	public static String polyAT_reverse_tag = "TT";

	public static String read_strand_tag = "RT";
	public static String flipped_tag="FT";

//	static String pAT_desc = "";
	public static String getHeader(){
		return "read_strand\tflipped\tdistA\tstA,endA\tdistT\tstT,endT";//+pAT_desc;
	}
public static String getInfo(SAMRecord sam) {
	Integer flip  = (Integer)sam.getAttribute(flipped_tag);
	String flipped;
	if(flip!=null){
		flipped = flip.intValue()==0 ? "unchanged": "flipped";
	}else{
		flipped = "null";
	}
	return sam.getAttribute(read_strand_tag)+"\t"+flipped+"\t"+sam.getAttribute(polyAT_forward_tag)+"\t"+sam.getAttribute(polyAT_reverse_tag);
}
	//
	
	static  int MODE = EdlibLibrary.EdlibAlignMode.EDLIB_MODE_HW;
	static  int TASK = EdlibLibrary.EdlibAlignTask.EDLIB_TASK_LOC;
//	static int maxsize_tolerance=1;
	static EdlibAlignConfig.ByValue config =EdlibLibrary.INSTANCE.edlibNewAlignConfig(-1, MODE, TASK, null, 0);
	
	// public PolyAT() throws IOException{
	//		}
	// need separate for dRNA and cDNA and cDNA single cell because of barcodes
	
	
	
	static String polyA, polyT;
	public static void set_bc_len(int bc_len_AT_1){
		char[] pA = new char[bc_len_AT_1];
		Arrays.fill(pA, 'A');
		polyA = new String(pA);
		char[] pT = new char[bc_len_AT_1];
		Arrays.fill(pT, 'T');
		polyT = new String(pT);
		bc_len_AT = bc_len_AT_1;
	} 
	
//	static int tolerance = 2;

	
	static int editDist(EdlibAlignResult resA_l ){
		return resA_l.editDistance+ Math.max(0,bc_len_AT-(resA_l.endLocations.getValue()-resA_l.startLocations.getValue()+1));
	}
	
	
	
	
	public static String check(String sequence_, int pos) {
		int dist = max_dist_to_end_AT;
		int offset = Math.max(0, pos-dist);
		String sequence = sequence_.substring(offset, Math.min(sequence_.length(), pos+dist));
		//String sequence = sam.getReadString();
		EdlibAlignResult.ByValue resT_l=		EdlibLibrary.INSTANCE.edlibAlign(polyT, polyT.length(), sequence, sequence.length(), config);
		int distT_l=	editDist(resT_l);
		EdlibAlignResult.ByValue resA_r=		EdlibLibrary.INSTANCE.edlibAlign(polyA,polyA.length(), sequence, sequence.length(), config);
		int distA_r=	editDist(resA_r);
		int mind = Math.min(distA_r, distT_l);
		//EdlibAlignResult.ByValue res = distA_r < distT_l ? resA_r : resT_l;
		String type = distA_r < distT_l ? "A" : "T";
		int st_p = distA_r < distT_l  ? resA_r.startLocations.getValue() : resT_l.startLocations.getValue();
		String str =  pos+","+(offset+st_p)+","+polyA.length()+","+type+","+mind;
		clear(resA_r); 
		clear(resT_l);
		return str;
	}
	
	
	 static void clear(ByValue resA_r) {
		EdlibLibrary.INSTANCE.edlibFreeAlignResult(resA_r);
		resA_r.clear();
		
	}
	/** 
	 * if we find polyA that means its forward_read
	 * if we find polyT that means its reverse_read
	 * does not modify sam record
	 * returns polyA edit distance
	 * forward means we assume it has polyA
	 * reverse forward is default so as not to in
	 * */
	public static int assign(SAMRecord sam, String sequence, boolean forward, boolean reverse_forward, int[] start_end) {
		int read_len = sequence.length();
		
		
		int len = max_dist_to_end_AT;
		if(read_len < len) return bc_len_AT;
		int offset = forward ? Math.max(0, read_len-len) : dist_to_ignore_AT;
		
		String sequence_1  = 
				forward ? sequence.substring(offset, read_len-dist_to_ignore_AT) : sequence.substring(dist_to_ignore_AT, Math.min(len, read_len)) ;
		String sequence_ ;
		if(reverse_forward){
			sequence_ = forward ?  SequenceUtil.reverseComplement(sequence_1) : sequence_1; // 
		}else{ // reverscomplement the not forward strand
			sequence_ = !forward ?  SequenceUtil.reverseComplement(sequence_1) : sequence_1; // 

		}
		
		String barc =reverse_forward ?  polyT : polyA; //forward ? polyA : polyT;
		EdlibAlignResult.ByValue resA_r=		EdlibLibrary.INSTANCE.edlibAlign(barc,bc_len_AT, sequence_, sequence_.length(), config);
		int distA_r=	editDist(resA_r);
//		int	 st_A = resA_r.startLocations.getValue()+offset;
//		int end_A = read_len -(resA_r.endLocations.getValue()+offset);
		int st_A, end_A;
		if(reverse_forward){
			if(forward){ // rev complemented end of sequence
				end_A = resA_r.startLocations.getValue()+dist_to_ignore_AT;
				st_A =  read_len  - (resA_r.endLocations.getValue()+dist_to_ignore_AT);
			}else{
				st_A = resA_r.startLocations.getValue()+dist_to_ignore_AT;
				end_A = read_len  - (resA_r.endLocations.getValue()+dist_to_ignore_AT);
			}
		}else{
			int offset1 = Math.max(0, read_len-len);
			if(forward){  // this is normal
				end_A = read_len -(resA_r.endLocations.getValue()+offset1);
				st_A =  resA_r.startLocations.getValue()+offset1;
			}else{ // rev complemented the start of the sequence
				st_A = read_len -(resA_r.endLocations.getValue()+offset1);
				end_A =  resA_r.startLocations.getValue()+offset1;
			}
		}
		
		String  pAT = distA_r+"\t"+st_A+","+(end_A);
		sam.setAttribute(forward ? PolyAT.polyAT_forward_tag: PolyAT.polyAT_reverse_tag, pAT);
		start_end[0] = st_A;
		start_end[1] = end_A;
		/*if(verbose && false){
//			System.err.println();
			System.err.println("no polyA or polyT "+distA_r);
			String seqRC =SequenceUtil.reverseComplement(sequence);
			EdlibAlignResult resT=		EdlibLibrary.INSTANCE.edlibAlign(polyT, bc_len_AT, sequence, sequence.length(), config);
			int distT=	editDist(resT);
			EdlibAlignResult resA=		EdlibLibrary.INSTANCE.edlibAlign(polyA, bc_len_AT, sequence, sequence.length(), config);
			int distA=	editDist(resA);
			EdlibAlignResult res = distA < distT ? resA : resT;
			int st = res.startLocations.getValue(); int end = read_len - res.endLocations.getValue();
			System.err.println(distA+","+distT+","+st+","+end);
		}*/
		clear(resA_r);
		return distA_r;
		
		
		
	}
	
}
