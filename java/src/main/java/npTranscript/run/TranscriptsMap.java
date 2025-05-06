package npTranscript.run;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.List;
import java.util.SortedMap;
import java.util.TreeMap;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;

// keeps a uniform list of transcripts
public class TranscriptsMap{
	
//	static  int MODE = EdlibLibrary.EdlibAlignMode.EDLIB_MODE_NW;
//	static  int TASK = EdlibLibrary.EdlibAlignTask.EDLIB_TASK_DISTANCE;
//	static int maxsize_tolerance=1;
//	static EdlibAlignConfig.ByValue config =EdlibLibrary.INSTANCE.edlibNewAlignConfig(-1, MODE, TASK, null, 0);
	
	
	//final boolean saveReadID;
	public TranscriptsMap(boolean addCounts){
		this.addcount = addCounts;
		if(addCounts) counts = new ArrayList<Integer>();

	}
	/** check freq is how often to check to make sure it matches */
	public void appendCounts(File inF, File outF, int check_freq ) throws IOException {
		BufferedReader in = new BufferedReader(new InputStreamReader(new GZIPInputStream((new FileInputStream(inF)))));
		PrintWriter out = new PrintWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(outF))));
		if(!addcount) throw new RuntimeException("!!");
		String str = "";
		for(int i=0; (str = in.readLine())!=null; i++){
			if(check_freq>0 && i%check_freq==0){
				String[] st = str.split("\t");
				if(!st[0].equals(this.cluster_ids.get(i))) throw new RuntimeException("!! not a match");
			}
			out.println(str+"\t"+this.counts.get(i));
		//	String[] st = str.split("\t");
			
		}
		counts.clear();
		this.addcount = false;
		in.close();
		out.close();
	}

	public TranscriptsMap(PrintWriter pw, File inputF, boolean addCounts) {
		//this.saveReadID = saveReadID;
		this.pw = pw;
		this.addcount = addCounts;
		if(addCounts) counts = new ArrayList<Integer>();
		if(inputF!=null){
			try{
			//	File inputF = new File(input);
				if(inputF.exists()){
					InputStream is = new FileInputStream(inputF);
					if(inputF.getName().endsWith(".gz" )) is = new GZIPInputStream(is);
			BufferedReader br = new BufferedReader(new InputStreamReader(is));
			String str = "";
			for(int i=0; (str=br.readLine())!=null; i++){
				String[] st = str.split("\t");
				cluster_id_map.put(st[0], i);
				this.cluster_ids.add(st[0]);
				if(addcount)this.counts.add(0);
				pw.println(addcount ? str.substring(0, str.lastIndexOf("\t")) : str);
			}
			br.close();
				}
			}catch(IOException exc){
				exc.printStackTrace();
			}
		}
	}
	
	PrintWriter pw = null;
	private List<Comparable>cluster_ids = new ArrayList<Comparable>();
	private List<Integer> counts = null;//new ArrayList<Integer>();
	
	private SortedMap<Comparable, Integer> cluster_id_map = new TreeMap<Comparable, Integer>();
	
	public Integer get(Comparable transcript){
		return this.cluster_id_map.get(transcript);
	}
	
	
	/*static int editDist(String umi, String umi1){
			EdlibAlignResult.ByValue resA_r=		EdlibLibrary.INSTANCE.edlibAlign(umi,umi.length(), umi1	, umi1.length(), config);
			int distA = resA_r.editDistance;
			clear(resA_r);
			return distA;
	}
	 static void clear(ByValue resA_r) {
			EdlibLibrary.INSTANCE.edlibFreeAlignResult(resA_r);
			resA_r.clear();
			
		}*/
	
	 static int edit_thresh = 2;
	 
	/*static class SuperUmi implements Comparable{
	
		public SuperUmi(String barcode, String umi){
			
			this.barcode = barcode;
			
				int umilen = umi.length();
				int len1 = Math.min(umilen, barcode.length()+12);
				this.umi = umilen > barcode.length() ? umi.substring(barcode.length(), len1) : null;
				this.polyT = umilen - len1;
	//		PolyAT.check(sequence_, pos)
		}
		String barcode;
		String umi;
		int polyT;
		public String toString(){
			return // useUMI ? 
					barcode+"_"+umi+"_"+polyT ;//: barcode;
		}
		
		@Override
		public boolean equals(Object o){
			return this.compareTo(o)==0;
		}
//boolean useUMI = false;
		@Override
		public int compareTo(Object o) {
		
			String barcode1 = ((SuperUmi)o).barcode;
			if(!barcode1.equals(barcode)){
				return barcode.compareTo(barcode1);
			}else{
//				if(!useUMI) return 0;else
				{
					String umi1 = ((SuperUmi)o).umi;
					int dist = editDist(umi1, umi);
					if(dist >edit_thresh){
						return umi.compareTo(umi1);
					}else{
						return 0;
					}
				}
			}
		}
		
	}*/
	
	 boolean addcount;
	//assumes is absent
	public Integer getNext(Comparable transcript,  boolean add){
		Integer trans_ind = cluster_ids.size();
		if(add){
			cluster_id_map.put(transcript,  trans_ind);
			this.cluster_ids.add(transcript);
			if(addcount) this.counts.add(0);
		
		}
		return trans_ind;
	}
	public void increment(int trans_ind){
		if(addcount)counts.set(trans_ind, counts.get(trans_ind)+1);
	}
	public Integer getAndAdd(String transcript,String[] line,
			int read_id_index, int chrom_index, int start_index, int end_index, int breaks_index
			) {

		
		Integer trans_ind = this.cluster_id_map.get(transcript);
		if(trans_ind==null){
			trans_ind = getNext(transcript,  true);
			if(pw!=null){
				System.err.println(transcript+"\t"+trans_ind);
				pw.println(transcript+"\t"+
						line[read_id_index]+"\t"+
						line[chrom_index]+"\t"+ line[start_index]+"\t"+line[end_index]+"\t"+line[breaks_index]
						);
			}
		}
		
		return trans_ind;
	}
	public void add(Comparable transcript, Integer trans_ind){
		cluster_id_map.put(transcript,  trans_ind);
		this.cluster_ids.add(transcript);
		if(addcount) this.counts.add(0);
	}
	public String[] getArray() {
		return cluster_ids.toArray(new String[0]);
	}
	public String[] getArray(int offset, int len1) {
		String[] res = new String[len1];
		for(int i=0; i<len1; i++){
			res[i] = this.cluster_ids.get(i+offset).toString();
		}
		return res;
	}
}