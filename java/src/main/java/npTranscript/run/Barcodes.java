package npTranscript.run;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Iterator;
import java.util.List;
import java.util.Set;
import java.util.SortedSet;
import java.util.TreeSet;
import java.util.zip.GZIPInputStream;

import com.sun.jna.ptr.IntByReference;

import edlib.EdlibAlignConfig;
import edlib.EdlibAlignResult;
import edlib.EdlibLibrary;
import edlib.EdlibAlignResult.ByValue;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.SequenceUtil;
import japsa.tools.seq.SequenceUtils;
import npTranscript.NW.PolyAT;
import npTranscript.cluster.TranscriptUtils;

public class Barcodes {
public static boolean verbose=false;	
	public static String barcode_forward_tag = "BF";
	public static String barcode_reverse_tag = "BR";
public static String umi_tag = "UM";
	public static String barcode_confidence_tag = "BI";
	//public  static String barcode_dist_tag="BD";
//	public  static String barcode_num_tag="BN";
//	public  static String barcode_left_tag="BL";
//	public  static String barcode_start_tag="CS";// note BS conflicts with another tag in SequenceUtiles
//	public  static String barcode_end_tag="BE";
//	public  static String barcode_forward_tag="BF";




	
	//String sequence = null;;
	

	public static void main(String[] args){
		try{
			String inf = "./long_read_UK_72hpi_adult_barcodes.tsv.gz";
			boolean forward = true;

			Barcodes barc = new Barcodes(inf);
			String bc = "AAACCCAAGCGTTAGG";
			String seq="ACGTTAATTCCGTTTTTTT";
			String suff="AAAAAAAAAAAAAAAAAAAAAAAAAA";
			String sequence = seq+bc+suff;
			SortedSet<Integer> mins = new TreeSet<Integer>();
			//mins.toArray(new Integer[0]);
			int minv = barc.calcMatches(mins, 0, sequence);
			int index=0;
			List<String>barcodes = new ArrayList<String>();
			//List<Boolean > forward = new ArrayList<Boolean> ();
			List<Integer> inds1 = new ArrayList<Integer>();
			barc.getAll(mins,  barcodes, inds1, forward);
			System.err.println(minv + " "+barcodes+" "+forward+" "+inds1);
		}catch(Exception exc){
			exc.printStackTrace();
		}
	}
	void getAll(Collection<Integer>inds, List<String> res, List<Integer>inds1, boolean forward_barcode){
		for(Iterator<Integer> it = inds.iterator(); it.hasNext();){
			int i = it.next();
			String barcode=this.barcodes_forward.get(i);//+"/"+this.barcodes_reverse.get(i);;
//					forward_barcode ? this.barcodes_forward[index].get(i): 
		//	boolean positive = this.rev_compl ? ((i%2) ==0) : true;
		//	if(!positive){
		//		barcode = SequenceUtil.reverseComplement(barcode);
		//	}
			inds1.add(i);
			//inds1.add(this.rev_compl ? Math.floorDiv(i,2) : i);
			res.add(barcode);
			//forward.add(positive);
		}
	}
 final List<String> barcodes_forward ;
 //final List<String> barcodes_reverse ;

 
 /*public void setSequence(String sequence){
	 this.sequence = sequence;
 }*/
 
static  int MODE = EdlibLibrary.EdlibAlignMode.EDLIB_MODE_HW;
static  int TASK = EdlibLibrary.EdlibAlignTask.EDLIB_TASK_LOC;
public static int tolerance_barcode =1;
public static int maxsize_tolerance_barcode=1;
public static int barcode_extent=150;
public static int barcode_ignore=0;

public final int bc_len_barcode;

static EdlibAlignConfig.ByValue config =EdlibLibrary.INSTANCE.edlibNewAlignConfig(-1, MODE, TASK, null, 0);
			
 
 public Barcodes(String infiles ) throws IOException{
	// this.rev_compl=incl_rev_compl;
	// barcodes_forward = new List[infiles.length];
	// barcodes_reverse= new List[infiles.length];

	// for(int i=0; i<infiles.length; i++){
		// if(infiles[i]!=null){
		 barcodes_forward = new ArrayList<String>();
		// barcodes_reverse = new ArrayList<String>();

		 File in = new File(infiles);
		 if(!in.exists()) throw new RuntimeException("does not exist "+in.getAbsolutePath());		if(in.exists()){
		InputStream fis = new FileInputStream(in);
		if(in.getName().endsWith(".gz")){
			fis = new GZIPInputStream(fis);
		}
		BufferedReader br = new BufferedReader(new InputStreamReader(fis));
		String st = "";
		while((st = br.readLine())!=null){
			String bc = st.split("-")[0];
			barcodes_forward.add(bc);
		//		barcodes_reverse.add(SequenceUtil.reverseComplement(bc));
//						TranscriptUtils.revC(bc));
		}
		}
		 //}
	// }
	 bc_len_barcode = barcodes_forward.get(0).length();
	}
	public int calcMatches(Set<Integer> mins, int index, String sequence){//, boolean forward_barcodes){
		//Set<Integer> mins = new TreeSet<Integer>();
		int minv=Integer.MAX_VALUE;
		List<String> barcodes_i =  barcodes_forward ;//: barcodes_reverse;
		int len = barcodes_i.size();
		//for(int j=0; j<sequence.length; j++){
		for(int i=0;i<len; i++){
				int dist_i = match(barcodes_i.get(i), sequence);
				if(dist_i <minv){
					mins.clear();
					minv = dist_i;
					mins.add(i);
				}else if(dist_i ==minv){
					mins.add(i);
				}
			
		}
		//}
		return minv;
	}
	
	private  static EdlibAlignResult.ByValue align(String sequence, String bc){
	 EdlibAlignResult.ByValue res=		EdlibLibrary.INSTANCE.edlibAlign(bc, bc.length(), sequence, sequence.length(), config);
	 
//		 EdlibAlignResult res=	EdlibLibrary.INSTANCE.edlibAlign(sequence, sequence.length(), bc, bc.length(), config);
		 return res;
	}
	private int match(String bc, String sequence) {
		if(sequence.length() < bc_len_barcode) return Integer.MAX_VALUE;
		 EdlibAlignResult.ByValue res=	align(sequence, bc);
		int out = res.editDistance+ Math.max(0,bc_len_barcode-(res.endLocations.getValue()-res.startLocations.getValue()+1));
		clear(res);
		return out;
	}
	public static int getPos(String bc, String sequence) {
	//	if(sequence.length() < bc_len_barcode) return Integer.MAX_VALUE;
		 EdlibAlignResult.ByValue res=	align(sequence, bc);
		int out = res.editDistance>3  ? -1 : res.endLocations.getValue();
	//	int out = res.editDistance+ Math.max(0,bc_len_barcode-(res.endLocations.getValue()-res.startLocations.getValue()+1));
		clear(res);
		return out;
	}
	
	 static void clear(ByValue resA_r) {
			EdlibLibrary.INSTANCE.edlibFreeAlignResult(resA_r);
			resA_r.clear();
			
		}
	
	 String fortag = "-5'_for";
	 String revtag = "+3'_revC";
			
/** finds barcodes within len bp of end and start, excluding the last and first chop bp 
 * 
 * assumes forward _read == reverse barcode at 3'end +polyA   and reverse_read == forward barcode at 5'end + polyT
 * 
 * sa is the original read uncorrected by minimap2

 * */
	 
	public int assign(SAMRecord sam, String sa, boolean forward_read, int[] startend) {
		Arrays.fill(startend, -1);
		boolean forward_barcode = !forward_read;
		int index = (Integer) sam.getAttribute(SequenceUtils.src_tag);
		List<String> barcodes_i = barcodes_forward;//forward_barcode  ? barcodes_forward : barcodes_reverse;
		if(barcodes_i.size()==0) return this.bc_len_barcode;
			int read_len = sa.length();
				// chop should be less then len
				if(barcode_extent < barcode_ignore) throw new RuntimeException("!!");
				String str1 = 
						forward_barcode ? sa.substring(barcode_ignore,Math.min(barcode_extent,read_len)) : sa.substring(Math.max(0, read_len-barcode_extent), read_len-barcode_ignore);
				int offset = forward_barcode ?  barcode_ignore : read_len-barcode_extent;
				
				String str = forward_barcode ? str1 : SequenceUtil.reverseComplement(str1);
				SortedSet<Integer> mins = new TreeSet<Integer>();
				int minv = calcMatches(mins, index,str );
				if(mins.size()==1 && minv <=tolerance_barcode){
					int inds1 = mins.first();
					String barc = this.barcodes_forward.get(inds1);
					 int st=0; int end=0;
					 EdlibAlignResult.ByValue res=	align(str, barc);
					// int dist=res.editDistance+ Math.max(0,bc_len_barcode-(res.endLocations.getValue()-res.startLocations.getValue()+1));
					
						
						//	 int en;
							
							
//						 sb.append(st);sb.append(",");sb.append(en);sb.append("\t");
						
					if(forward_barcode){ //normal
						 st = res.startLocations.getValue()+barcode_ignore;
						 end = read_len - (res.endLocations.getValue()+barcode_ignore);
					}else{//rev complemented end of sequence
						end = res.startLocations.getValue()+barcode_ignore;
						st = read_len - (res.endLocations.getValue()+barcode_ignore);
					}
					 String sb = barc+"\t"+inds1+"\t"+minv+"\t"+st+","+end+"\t"+(forward_barcode ? fortag: revtag);
					 startend[0] = st;
					 startend[1] = end;
						 
						 
							 clear(res);
						sam.setAttribute(forward_read ? Barcodes.barcode_forward_tag : Barcodes.barcode_reverse_tag, sb);
						return minv;
				}else{
					if( verbose) System.err.println("not found "+minv+" "+mins.size());
				}
				return minv;
	}
	public static String getHeader(){
		return "barcode\tbarcode_index\tbarcode_dist\tbarcode_pos\tbarcode_type\tconfidence\tUMI";
	}
	static String no_res = "null\tNA\tNA\tnull\tnone";
	static String both_res = "null\tNA\tNA\tnull\tboth";
public static String getInfo(SAMRecord sam) {
	String  forward = (String) sam.getAttribute(barcode_forward_tag);
	String rev = (String) sam.getAttribute(barcode_reverse_tag);
	String comb = no_res;
	String conf = (String) sam.getAttribute(barcode_confidence_tag);;
	String umi = (String) sam.getAttribute(umi_tag);;
	if(forward==null && rev==null){
		
	   comb = no_res;	
	}else if(forward!=null && rev!=null){
		String[] str1 = forward.split("\t");
		String[] str2 = rev.split("\t");
		int f = Integer.parseInt(str1[2]);
		int r =  Integer.parseInt(str2[2]);
		if(f<r){
			comb=forward;
			conf=conf+",both,"+f+","+r;
		}else if(r<f){
			comb = rev;
			conf=conf+",both,"+f+","+r;
		}else {
			comb = both_res+","+f+","+r;
		}
	}else if(forward!=null){
		comb = forward;
	}else if(rev!=null){
		comb = rev;
	}
	return comb +"\t"+conf+"\t"+umi;
			
}
	
}
