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
import edlib.EdlibLibrary;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.SequenceUtil;
import npTranscript.cluster.TranscriptUtils;

public class PolyAT {
	
	public static String polyAT_tag = "AT";
	public static String read_strand_tag = "RT";

	static  int MODE = EdlibLibrary.EdlibAlignMode.EDLIB_MODE_HW;
	static  int TASK = EdlibLibrary.EdlibAlignTask.EDLIB_TASK_LOC;
	static int tolerance = 2;
	static int maxsize_tolerance=1;
	static EdlibAlignConfig.ByValue config =EdlibLibrary.INSTANCE.edlibNewAlignConfig(-1, MODE, TASK, null, 0);
	
	// public PolyAT() throws IOException{
	//		}
	
	static final int bc_len = 20;
	static String polyA, polyT;
	static{
		char[] pA = new char[bc_len];
		Arrays.fill(pA, 'A');
		polyA = new String(pA);
		char[] pT = new char[bc_len];
		Arrays.fill(pT, 'T');
		polyT = new String(pT);
	} 
	static int edit_thresh=2;
	static int max_dist_to_end=100;
	
	static int editDist(EdlibAlignResult resA_l ){
		return resA_l.editDistance+ Math.max(0,bc_len-(resA_l.endLocations.getValue()-resA_l.startLocations.getValue()+1));
	}
	/** 
	 * if we find polyA that means its forward_read
	 * if we find polyT that means its reverse_read
	 * does not modify sam record
	 * returns if this read strand is positive (ie if it has a polyA 
	 * */
	public Boolean assign(SAMRecord sam) {
		String sequence = sam.getReadString();
		int read_len = sequence.length();
		boolean neg = sam.getReadNegativeStrandFlag();
		 if(neg){ // this means the sequence has been rev complemented, so this corrects back to the original read from fastq
			 sequence = SequenceUtil.reverseComplement(sequence);
		 }
		int len = max_dist_to_end;
		if(read_len < len) return null;
		int end_offset = Math.max(0, read_len-len);
		int start_offset = 0;
		String sequenceL = sequence.substring(0,Math.min(len,read_len));
		String sequenceR = sequence.substring(end_offset, read_len);

		//EdlibAlignResult resA_l=		EdlibLibrary.INSTANCE.edlibAlign(polyA,bc_len, sequenceL, sequenceL.length(), config);
		//int distA_l=	editDist(resA_l);
		//EdlibAlignResult resT_r=		EdlibLibrary.INSTANCE.edlibAlign(polyT, bc_len, sequenceR, sequenceR.length(), config);
				//int distT_r=	editDist(resT_r);
		
		EdlibAlignResult resT_l=		EdlibLibrary.INSTANCE.edlibAlign(polyT, bc_len, sequenceL, sequenceL.length(), config);
		int distT_l=	editDist(resT_l);
		EdlibAlignResult resA_r=		EdlibLibrary.INSTANCE.edlibAlign(polyA,bc_len, sequenceR, sequenceR.length(), config);
		int distA_r=	editDist(resA_r);

		
		
		//int distA = Math.min(distA_l, distA_r);
		//int distT = Math.min(distT_l, distT_r);

		//need the right distance less than left_distance
		if(distA_r<=edit_thresh && distA_r<distT_l) { // && distA_r < distA_l){
			// this means the original sequence is forward strand
			int	 st_A = resA_r.startLocations.getValue()+end_offset;
			int end_A = read_len -(resA_r.endLocations.getValue()+end_offset);
			String  pAT = "A,"+st_A+","+(end_A)+","+distA_r;
			 sam.setAttribute(read_strand_tag, "+");
			 sam.setAttribute(PolyAT.polyAT_tag, pAT);
/*			if(end_A > max_dist_to_end){
				return null;
			}*/
				
			 System.err.println(distA_r+" "+distT_l+" "+pAT);
			 return true;
		}else if(distT_l <= edit_thresh  && distT_l < distA_r){// && distT_l < distT_r){
			// this means the original sequence is reverse strand
			
			int	 st_T = resT_l.startLocations.getValue() ;
			int end_T = read_len -resT_l.endLocations.getValue();
			String pAT = "T,"+st_T+","+(end_T)+","+distT_l;
			 sam.setAttribute(PolyAT.read_strand_tag, "-");
			 sam.setAttribute(PolyAT.polyAT_tag, pAT);
			 //should be at 5'end
			// if(Math.min(end_A, st_A)> max_dist_to_end)  return null;
			/* if(st_T > max_dist_to_end) {
				 return null;
			 }*/
			 System.err.println(distA_r+" "+distT_l+" "+pAT);
			 return false;

		}else{
		//	System.err.println(distA_r+" "+distT_l);
		//	System.err.println("no polyA or polyT");

		}
		return null;
		
		
		
	}
}
