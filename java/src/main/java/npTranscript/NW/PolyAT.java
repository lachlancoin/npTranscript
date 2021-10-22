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
import japsa.seq.Sequence;
import npTranscript.cluster.TranscriptUtils;

public class PolyAT {
	private static boolean verbose=false;
	private static  int bc_len_AT = 6; //perhaps should be 15. seems shorter for RNA
	public static int edit_thresh_AT=1; // how many mismatches
	public static int max_dist_to_end_AT=20;  //how far from ends to search  150 for cDNA 100 for dRNA
	
	public static String polyAT_tag = "AT";
	public static String read_strand_tag = "RT";
	public static String flipped_tag="FT";

//	static String pAT_desc = "";
	public static String getHeader(){
		return "read_strand\tA/T,st,end,dist";//+pAT_desc;
	}
public static String getInfo(SAMRecord sam) {
	return sam.getAttribute(read_strand_tag)+"\t"+sam.getAttribute(polyAT_tag);
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
		EdlibAlignResult resT_l=		EdlibLibrary.INSTANCE.edlibAlign(polyT, polyT.length(), sequence, sequence.length(), config);
		int distT_l=	editDist(resT_l);
		EdlibAlignResult resA_r=		EdlibLibrary.INSTANCE.edlibAlign(polyA,polyA.length(), sequence, sequence.length(), config);
		int distA_r=	editDist(resA_r);
		int mind = Math.min(distA_r, distT_l);
		EdlibAlignResult res = distA_r < distT_l ? resA_r : resT_l;
		String type = distA_r < distT_l ? "A" : "T";
		return pos+","+(offset+res.startLocations.getValue())+","+polyA.length()+","+type+","+mind;
		//return distA_r+","+distT_l;
/*		if(mind<=edit_thresh_AT){
			
		}else{
			return null;
		}*/
	}
	
	
	/** 
	 * if we find polyA that means its forward_read
	 * if we find polyT that means its reverse_read
	 * does not modify sam record
	 * returns if this read strand is positive (ie if it has a polyA 
	 * */
	public static Boolean assign(SAMRecord sam) {
		String sequence = sam.getReadString();
		int read_len = sequence.length();
		boolean neg = sam.getReadNegativeStrandFlag();
		 if(neg){ // this means the sequence has been rev complemented, so this corrects back to the original read from fastq
			 sequence = SequenceUtil.reverseComplement(sequence);
		 }
		int len = max_dist_to_end_AT;
		if(read_len < len) return null;
		int end_offset = Math.max(0, read_len-len);
		int start_offset = 0;
		String sequenceL = sequence.substring(0,Math.min(len,read_len));
		String sequenceR = sequence.substring(end_offset, read_len);

		//EdlibAlignResult resA_l=		EdlibLibrary.INSTANCE.edlibAlign(polyA,bc_len, sequenceL, sequenceL.length(), config);
		//int distA_l=	editDist(resA_l);
		//EdlibAlignResult resT_r=		EdlibLibrary.INSTANCE.edlibAlign(polyT, bc_len, sequenceR, sequenceR.length(), config);
				//int distT_r=	editDist(resT_r);
		
		EdlibAlignResult resT_l=		EdlibLibrary.INSTANCE.edlibAlign(polyT, bc_len_AT, sequenceL, sequenceL.length(), config);
		int distT_l=	editDist(resT_l);
		EdlibAlignResult resA_r=		EdlibLibrary.INSTANCE.edlibAlign(polyA,bc_len_AT, sequenceR, sequenceR.length(), config);
		int distA_r=	editDist(resA_r);

		
		
		//int distA = Math.min(distA_l, distA_r);
		//int distT = Math.min(distT_l, distT_r);

		//need the right distance less than left_distance
		if(distA_r<=edit_thresh_AT && distA_r<distT_l) { // && distA_r < distA_l){
			// this means the original sequence is forward strand
			int	 st_A = resA_r.startLocations.getValue()+end_offset;
			int end_A = read_len -(resA_r.endLocations.getValue()+end_offset);
			String  pAT = "A,"+st_A+","+(end_A)+","+distA_r;
			 sam.setAttribute(read_strand_tag, "+");
			 sam.setAttribute(PolyAT.polyAT_tag, pAT);
/*			if(end_A > max_dist_to_end){
				return null;
			}*/
				
		//	 System.err.println(distA_r+" "+distT_l+" "+pAT);
			 return true;
		}else if(distT_l <= edit_thresh_AT  && distT_l < distA_r){// && distT_l < distT_r){
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
			// System.err.println(distA_r+" "+distT_l+" "+pAT);
			 return false;

		}else if(verbose){
//			System.err.println();
			System.err.println("no polyA or polyT "+distA_r+" "+distT_l);
			String seqRC =SequenceUtil.reverseComplement(sequence);
			EdlibAlignResult resT=		EdlibLibrary.INSTANCE.edlibAlign(polyT, bc_len_AT, sequence, sequence.length(), config);
			int distT=	editDist(resT);
			EdlibAlignResult resA=		EdlibLibrary.INSTANCE.edlibAlign(polyA, bc_len_AT, sequence, sequence.length(), config);
			int distA=	editDist(resA);
			EdlibAlignResult res = distA < distT ? resA : resT;
			int st = res.startLocations.getValue(); int end = read_len - res.endLocations.getValue();
			System.err.println(distA+","+distT+","+st+","+end);
		}
		return null;
		
		
		
	}
	
}
