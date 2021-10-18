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

import edlib.EdlibAlignConfig;
import edlib.EdlibAlignResult;
import edlib.EdlibLibrary;
import htsjdk.samtools.SAMRecord;
import japsa.tools.seq.SequenceUtils;

public class Barcodes {
	
	public static String barcode_tag = "BC";
	public  static String barcode_dist_tag="BD";
	public  static String barcode_num_tag="BN";
	public  static String barcode_left_tag="BL";

	
	String sequence = null;;
	

	public static void main(String[] args){
		try{
			String inf = "./long_read_UK_72hpi_adult_barcodes.tsv.gz";
			Barcodes barc = new Barcodes(new String[] {inf});
			String bc = "AAACCCAAGCGTTAGG";
			String seq="ACGTTAATTCCGTTTTTTT";
			String suff="AAAAAAAAAAAAAAAAAAAAAAAAAA";
			barc.setSequence(seq+bc+suff);
			SortedSet<Integer> mins = new TreeSet<Integer>();
			//mins.toArray(new Integer[0]);
			int minv = barc.calcMatches(mins);
			int index=0;
			System.err.println(minv + " "+barc.getAll(mins, index));
		}catch(Exception exc){
			exc.printStackTrace();
		}
	}
	List<String> getAll(Collection<Integer>inds, int index){
		List<String> res = new ArrayList<String>();
		for(Iterator<Integer> it = inds.iterator(); it.hasNext();){
			res.add(this.barcodes[index].get(it.next()));
		}
		return res;
	}
 final List<String>[] barcodes ;
 public void setSequence(String sequence){
	 this.sequence = sequence;
 }
 
static  int MODE = EdlibLibrary.EdlibAlignMode.EDLIB_MODE_HW;
static  int TASK = EdlibLibrary.EdlibAlignTask.EDLIB_TASK_LOC;
static int tolerance = 2;
static int maxsize_tolerance=1;
static EdlibAlignConfig.ByValue config =EdlibLibrary.INSTANCE.edlibNewAlignConfig(-1, MODE, TASK, null, 0);
			
 
 
 public Barcodes(String[] infiles ) throws IOException{
	 barcodes = new List[infiles.length];
	 for(int i=0; i<infiles.length; i++){
		 barcodes[i] = new ArrayList<String>();
		 File in = new File(infiles[i]);
		if(in.exists()){
		InputStream fis = new FileInputStream(in);
		if(in.getName().endsWith(".gz")){
			fis = new GZIPInputStream(fis);
		}
		BufferedReader br = new BufferedReader(new InputStreamReader(fis));
		String st = "";
		while((st = br.readLine())!=null){
			barcodes[i].add(st.split("-")[0]);
		}
		}
	 }
	}
	public int calcMatches(Set<Integer> mins, int index){
		//Set<Integer> mins = new TreeSet<Integer>();
		int minv=Integer.MAX_VALUE;
		List<String> barcodes_i = barcodes[index];
		int len = barcodes_i.size();
		for(int i=0;i<len; i++){
			int dist_i = match(barcodes_i.get(i));
			if(dist_i <minv){
				mins.clear();
				minv = dist_i;
				mins.add(i);
			}else if(dist_i ==minv){
				mins.add(i);
			}
		}
		return minv;
	}
	private  EdlibAlignResult align(String sequence, String bc){
	 EdlibAlignResult res=		EdlibLibrary.INSTANCE.edlibAlign(bc, bc.length(), sequence, sequence.length(), config);
//		 EdlibAlignResult res=	EdlibLibrary.INSTANCE.edlibAlign(sequence, sequence.length(), bc, bc.length(), config);
		 return res;
	}
	private int match(String bc) {
		// TODO Auto-generated method stub
		 // StringByReference target = new StringByReference(bc);
		 EdlibAlignResult res=	align(sequence, bc);
		 int st = res.startLocations.getValue();
		 int end = res.endLocations.getValue();
		 int len = end-st+1;
		
		 int dist=res.editDistance;
 if(len < bc.length()){
	 int diff = bc.length()-len;
			dist += diff ;
		 }
		/* System.err.println("start: "+st);
		 System.err.println("end: "+ end);
		 System.err.println("length:"+  len);
		 System.err.println("dist:"+  dist);*/
		return dist;
	}
	
	
/** finds barcodes within len bp of end and start, excluding the last and first chop bp */
	public void assign(SAMRecord sam,int len, int chop) {
		int index = (Integer) sam.getAttribute(SequenceUtils.src_tag);
		List<String> barcodes_i = barcodes[index];
		if(barcodes_i.size()>0){
			int read_len = sam.getReadLength();
			if(read_len<=2*len){
				this.setSequence(sam.getReadString().substring(chop, read_len-chop));
				SortedSet<Integer> mins = new TreeSet<Integer>();
				int minv = calcMatches(mins, index);
				if(minv <=tolerance && mins.size()<=maxsize_tolerance){
						sam.setAttribute(barcode_tag, mins.first());
						sam.setAttribute(barcode_dist_tag, minv);
						sam.setAttribute(barcode_num_tag, mins.size());
				}
			}else{
				this.setSequence(sam.getReadString().substring(chop,len));
				SortedSet<Integer> mins_l = new TreeSet<Integer>();
				int minv_l = calcMatches(mins_l, index);
				
				this.setSequence(sam.getReadString().substring(read_len-len, read_len-chop));
				SortedSet<Integer> mins_r = new TreeSet<Integer>();
				int minv_r = calcMatches(mins_r, index);
				int minv = Math.min(minv_l, minv_r);
				SortedSet<Integer>  mins = minv_l < minv_r ? mins_l : mins_r;
				boolean left_tag = minv_l < minv_r ? true: false;
				if(minv <=tolerance && mins.size()<=maxsize_tolerance){
					sam.setAttribute(barcode_tag, mins.first());
					sam.setAttribute(barcode_dist_tag, minv);
					sam.setAttribute(barcode_num_tag, mins.size());
					sam.setAttribute(barcode_left_tag, left_tag);
				}
			}
		}
	}
	
}
