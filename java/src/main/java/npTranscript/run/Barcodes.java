package npTranscript.run;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
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

			Barcodes barc = new Barcodes(new String[] {inf});
			String bc = "AAACCCAAGCGTTAGG";
			String seq="ACGTTAATTCCGTTTTTTT";
			String suff="AAAAAAAAAAAAAAAAAAAAAAAAAA";
			String sequence = seq+bc+suff;
			SortedSet<Integer> mins = new TreeSet<Integer>();
			//mins.toArray(new Integer[0]);
			int minv = barc.calcMatches(mins, 0, sequence, forward);
			int index=0;
			List<String>barcodes = new ArrayList<String>();
			//List<Boolean > forward = new ArrayList<Boolean> ();
			List<Integer> inds1 = new ArrayList<Integer>();
			barc.getAll(mins, index, barcodes, inds1, forward);
			System.err.println(minv + " "+barcodes+" "+forward+" "+inds1);
		}catch(Exception exc){
			exc.printStackTrace();
		}
	}
	void getAll(Collection<Integer>inds, int index, List<String> res, List<Integer>inds1, boolean forward_barcode){
		for(Iterator<Integer> it = inds.iterator(); it.hasNext();){
			int i = it.next();
			String barcode=forward_barcode ? this.barcodes_forward[index].get(i): this.barcodes_reverse[index].get(i);
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
 final List<String>[] barcodes_forward ;
 final List<String>[] barcodes_reverse ;

 
 /*public void setSequence(String sequence){
	 this.sequence = sequence;
 }*/
 
static  int MODE = EdlibLibrary.EdlibAlignMode.EDLIB_MODE_HW;
static  int TASK = EdlibLibrary.EdlibAlignTask.EDLIB_TASK_LOC;
public static int tolerance_barcode = 2;
public static int maxsize_tolerance_barcode=1;
public static int barcode_extent=150;
public static int barcode_ignore=0;

public final int bc_len_barcode;

static EdlibAlignConfig.ByValue config =EdlibLibrary.INSTANCE.edlibNewAlignConfig(-1, MODE, TASK, null, 0);
			
 
 public Barcodes(String[] infiles ) throws IOException{
	// this.rev_compl=incl_rev_compl;
	 barcodes_forward = new List[infiles.length];
	 barcodes_reverse= new List[infiles.length];

	 for(int i=0; i<infiles.length; i++){
		 if(infiles[i]!=null){
		 barcodes_forward[i] = new ArrayList<String>();
		 barcodes_reverse[i] = new ArrayList<String>();

		 File in = new File(infiles[i]);
		 if(!in.exists()) throw new RuntimeException("does not exist "+in.getAbsolutePath());		if(in.exists()){
		InputStream fis = new FileInputStream(in);
		if(in.getName().endsWith(".gz")){
			fis = new GZIPInputStream(fis);
		}
		BufferedReader br = new BufferedReader(new InputStreamReader(fis));
		String st = "";
		while((st = br.readLine())!=null){
			String bc = st.split("-")[0];
			barcodes_forward[i].add(bc);
				barcodes_reverse[i].add(SequenceUtil.reverseComplement(bc));
//						TranscriptUtils.revC(bc));
		}
		}
		 }
	 }
	 bc_len_barcode = barcodes_forward[0].get(0).length();
	}
	public int calcMatches(Set<Integer> mins, int index, String sequence, boolean forward_barcodes){
		//Set<Integer> mins = new TreeSet<Integer>();
		int minv=Integer.MAX_VALUE;
		List<String> barcodes_i = forward_barcodes ? barcodes_forward[index] : barcodes_reverse[index];
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
	
	 static void clear(ByValue resA_r) {
			EdlibLibrary.INSTANCE.edlibFreeAlignResult(resA_r);
			resA_r.clear();
			
		}
	
	
			
/** finds barcodes within len bp of end and start, excluding the last and first chop bp 
 * 
 * assumes forward _read == reverse barcode at 3'end  and reverse_read == forward barcode at 5'end
 * sa is the original read uncorrected by minimap2
 * */
	public int assign(SAMRecord sam, String sa, boolean forward_read) {
	//	boolean neg_align=sam.getReadNegativeStrandFlag();
    //	boolean forward_align = !neg_align; // whether read is forward, we need this so we know if we have to rev compl to get original orientation
		
		/*String read_strand_tag = (String) sam.getAttribute(PolyAT.read_strand_tag);
		if(read_strand_tag==null){
			return false ;
		}
		boolean forward_read = read_strand_tag.charAt(0)=='+';*/
		boolean forward_barcode = !forward_read;
		

		
		int index = (Integer) sam.getAttribute(SequenceUtils.src_tag);
		List<String> barcodes_i = forward_barcode  ? barcodes_forward[index] : barcodes_reverse[index];
		if(barcodes_i.size()==0) return this.bc_len_barcode;
			
			int read_len = sa.length();
				// chop should be less then len
				if(barcode_extent < barcode_ignore) throw new RuntimeException("!!");
				String sequenceL = sa.substring(barcode_ignore,Math.min(barcode_extent,read_len));
			//	System.err.println(read_len+" "+sequenceL);
				//if(true) System.exit(0);;
				//System.err.println(read_len+" "+chop);
				String sequenceR = sa.substring(Math.max(0, read_len-barcode_extent), read_len-barcode_ignore);
				String str = forward_read ?  sequenceR  :  sequenceL;
				int offset = forward_read ? read_len-barcode_extent : barcode_ignore;
		//		String[] desc = forward_read ? new String[] {"3'"} : new String[] {"5'"} ;

				SortedSet<Integer> mins = new TreeSet<Integer>();
				int minv = calcMatches(mins, index,str, forward_barcode );
				if(mins.size()<=maxsize_tolerance_barcode){
					List<String>barcodes1 = new ArrayList<String>();
				//	List<Boolean > forward_l = new ArrayList<Boolean> (); // whether barcode is forward
					List<Integer> inds1 = new ArrayList<Integer>();
					getAll(mins, index, barcodes1, inds1, forward_barcode);
				//	String left = "NA";
					
					//sam.setAttribute(barcode_num_tag, mins.size());
					
					 int st=0; int end=0;
					if(barcodes1.size()==1 && minv <=tolerance_barcode){
						StringBuffer sb = new StringBuffer();
						
						sb.append(mins.size()==1 ? barcodes1.get(0) : barcodes1.toString());sb.append("\t");
						sb.append( mins.size()==1 ? inds1.get(0) :  inds1.toString());sb.append("\t");
						sb.append(minv);sb.append("\t");
						String barc = barcodes1.get(0);
						 int min_dist=Integer.MAX_VALUE;
						// int min_ind=0;
						
					//	for(int j=0; j<str.length; j++){
							EdlibAlignResult.ByValue res=	align(str, barc);
							 int dist=res.editDistance+ Math.max(0,bc_len_barcode-(res.endLocations.getValue()-res.startLocations.getValue()+1));
							 if(dist <min_dist){
								 st = res.startLocations.getValue()+offset;
								 end = res.endLocations.getValue()+offset;
							//	 min_ind=j;
								 min_dist = dist;
							 }
							 clear(res);
//							 res.clear();
					//	}
						//left=min_ind==0? "5'":"3'";
						 //sb.append(desc[min_ind]);sb.append(",");
						 sb.append(st);sb.append(",");sb.append(read_len-end);sb.append("\t");
						//sb.append(forward_barcode ? "forward": "revC");
						//Object pAT = sam.getAttribute(PolyAT.polyAT_tag);
						
				//if(verbose)	System.err.println((forward_align? "align+" : "align-")+" "+(forward_read? "read+" : "read-")+" "+minv + " "+barcodes1+" "
					//	 +(forward_barcode ? "forward_barcode": "rev_barcode")+" "+inds1+" "+desc[min_ind]+ " "+st+" "+(read_len-end)+(forward_read ? " polyA:" : " polyT:")+pAT);
						if(forward_read) sb.append("+3'_revC"); else sb.append("-5'_for");
						
						sam.setAttribute(forward_read ? Barcodes.barcode_forward_tag : Barcodes.barcode_reverse_tag, sb.toString());
						return minv;
					}
				}else{
					if( verbose) System.err.println("not found "+minv+" "+mins.size());
				}
				return minv;
	}
	public static String getHeader(){
		return "barcode\tbarcode_index\tbarcode_dist\tbarcode_pos\tbarcode_type\tconfidence";
	}
	static String no_res = "null\tNA\tNA\tnull\tnone";
	static String both_res = "null\tNA\tNA\tnull\tboth";
public static String getInfo(SAMRecord sam) {
	String  forward = (String) sam.getAttribute(barcode_forward_tag);
	String rev = (String) sam.getAttribute(barcode_reverse_tag);
	String comb = no_res;
	String conf = (String) sam.getAttribute(barcode_confidence_tag);;
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
	return comb +"\t"+conf;
			
}
	
}
