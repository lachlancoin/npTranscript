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
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.SequenceUtil;
import japsa.tools.seq.SequenceUtils;
import npTranscript.NW.PolyAT;
import npTranscript.cluster.TranscriptUtils;

public class Barcodes {
	
	public static String barcode_tag = "BC";
	public static String barcode_index_tag = "BI";
	public  static String barcode_dist_tag="BD";
	public  static String barcode_num_tag="BN";
	public  static String barcode_left_tag="BL";
	public  static String barcode_start_tag="CS";// note BS conflicts with another tag in SequenceUtiles
	public  static String barcode_end_tag="BE";
	public  static String barcode_forward_tag="BF";




	
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
			int minv = barc.calcMatches(mins, 0, new String[] {sequence}, forward);
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
static int tolerance = 2;
static int maxsize_tolerance=1;
static EdlibAlignConfig.ByValue config =EdlibLibrary.INSTANCE.edlibNewAlignConfig(-1, MODE, TASK, null, 0);
			
/** assumes bc_len is constant */
 final int bc_len;
 
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
				barcodes_reverse[i].add(TranscriptUtils.revC(bc));
		}
		}
		 }
	 }
	 bc_len = barcodes_forward[0].get(0).length();
	}
	public int calcMatches(Set<Integer> mins, int index, String[] sequence, boolean forward_barcodes){
		//Set<Integer> mins = new TreeSet<Integer>();
		int minv=Integer.MAX_VALUE;
		List<String> barcodes_i = forward_barcodes ? barcodes_forward[index] : barcodes_reverse[index];
		int len = barcodes_i.size();
		for(int j=0; j<sequence.length; j++){
		for(int i=0;i<len; i++){
				int dist_i = match(barcodes_i.get(i), sequence[j]);
				if(dist_i <minv){
					mins.clear();
					minv = dist_i;
					mins.add(i);
				}else if(dist_i ==minv){
					mins.add(i);
				}
			
		}
		}
		return minv;
	}
	
	private  static EdlibAlignResult align(String sequence, String bc){
	 EdlibAlignResult res=		EdlibLibrary.INSTANCE.edlibAlign(bc, bc.length(), sequence, sequence.length(), config);
//		 EdlibAlignResult res=	EdlibLibrary.INSTANCE.edlibAlign(sequence, sequence.length(), bc, bc.length(), config);
		 return res;
	}
	private int match(String bc, String sequence) {
		if(sequence.length() < bc_len) return Integer.MAX_VALUE;
		 EdlibAlignResult res=	align(sequence, bc);
		 return res.editDistance+ Math.max(0,bc_len-(res.endLocations.getValue()-res.startLocations.getValue()+1));
	}
	
	
			
/** finds barcodes within len bp of end and start, excluding the last and first chop bp */
	public void assign(SAMRecord sam,int len, int chop) {
		boolean neg_align=sam.getReadNegativeStrandFlag();
		
		
		
	boolean forward_align = !neg_align; // whether read is forward
		
		String read_strand_tag = (String) sam.getAttribute(PolyAT.read_strand_tag);
		if(read_strand_tag==null){
			return ;
		}
		boolean forward_read = read_strand_tag.charAt(0)=='+';
		boolean forward_barcode = !forward_read;
		

		
		int index = (Integer) sam.getAttribute(SequenceUtils.src_tag);
		List<String> barcodes_i = forward_barcode  ? barcodes_forward[index] : barcodes_reverse[index];
		if(barcodes_i.size()==0) return;
			String sa = sam.getReadString();
			if(neg_align){
				// this converts read back to original orientation
				sa = SequenceUtil.reverseComplement(sa);
			}
			int read_len = sa.length();
			{
				// chop should be less then len
				if(len < chop) throw new RuntimeException("!!");
				String sequenceL = sa.substring(chop,Math.min(len,read_len));
			//	System.err.println(read_len+" "+sequenceL);
				//if(true) System.exit(0);;
				//System.err.println(read_len+" "+chop);
				String sequenceR = sa.substring(Math.max(0, read_len-len), read_len-chop);
				String[] str = forward_read ?  new String[] {sequenceR}  :  new String[] {sequenceL};
				int[] offset = forward_read ? new int[] {read_len-len} : new int[] {chop};
				String[] desc = forward_read ? new String[] {"3'"} : new String[] {"5'"} ;

				SortedSet<Integer> mins = new TreeSet<Integer>();
				int minv = calcMatches(mins, index,str, forward_barcode );
				if(mins.size()<=maxsize_tolerance){
					List<String>barcodes1 = new ArrayList<String>();
				//	List<Boolean > forward_l = new ArrayList<Boolean> (); // whether barcode is forward
					List<Integer> inds1 = new ArrayList<Integer>();
					getAll(mins, index, barcodes1, inds1, forward_barcode);
				//	String left = "NA";
					
					sam.setAttribute(barcode_tag, barcodes1.toString());
					sam.setAttribute(barcode_index_tag, inds1.toString());
					sam.setAttribute(barcode_dist_tag, minv);
					sam.setAttribute(barcode_num_tag, mins.size());
					
					 int st=0; int end=0;
					if(barcodes1.size()==1 && minv <=tolerance){
						String barc = barcodes1.get(0);
						 int min_dist=Integer.MAX_VALUE;
						 int min_ind=0;
						
						for(int j=0; j<str.length; j++){
							EdlibAlignResult res=	align(str[j], barc);
							 int dist=res.editDistance+ Math.max(0,bc_len-(res.endLocations.getValue()-res.startLocations.getValue()+1));
							 if(dist <min_dist){
								 st = res.startLocations.getValue()+offset[j];
								 end = res.endLocations.getValue()+offset[j];
								 min_ind=j;
								 min_dist = dist;
							 }
						}
						//left=min_ind==0? "5'":"3'";
						 sam.setAttribute(barcode_left_tag,desc[min_ind]);
						 sam.setAttribute(barcode_start_tag, st);
						 sam.setAttribute(barcode_end_tag, read_len-end);
						 sam.setAttribute(barcode_forward_tag,forward_barcode ? "forward": "revC");
						Object pAT = sam.getAttribute(PolyAT.polyAT_tag);
						
					System.err.println((forward_align? "align+" : "align-")+" "+(forward_read? "read+" : "read-")+" "+minv + " "+barcodes1+" "
						 +(forward_barcode ? "forward_barcode": "rev_barcode")+" "+inds1+" "+desc[min_ind]+ " "+st+" "+(read_len-end)+(forward_read ? " polyA:" : " polyT:")+pAT);

					}
				}else{
					System.err.println("not found "+minv+" "+mins.size());
				}
			}
	}
	public static String getHeader(){
		return "barcode\tbarcode_index\tbarcode_dist\tbarcode_num\tbarcode_left\tbarcode_forward\tdist_from_start\tdist_from_end";
	}
public static String getInfo(SAMRecord sam) {
	return sam.getAttribute(barcode_tag)+"\t"+sam.getAttribute(barcode_index_tag)+"\t"+sam.getAttribute(barcode_dist_tag)
	+"\t"+sam.getAttribute(barcode_num_tag)+"\t"+sam.getAttribute(barcode_left_tag)+"\t"+sam.getAttribute(barcode_forward_tag)+"\t"+sam.getAttribute(barcode_start_tag)+"\t"+sam.getAttribute(barcode_end_tag);
}
	
}
