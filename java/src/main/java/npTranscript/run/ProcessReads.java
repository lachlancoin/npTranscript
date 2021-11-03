package npTranscript.run;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.OutputStream;
import java.io.OutputStreamWriter;
import java.io.PrintWriter;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import ch.systemsx.cisd.hdf5.IHDF5Writer;
import htsjdk.samtools.util.SequenceUtil;
import npTranscript.cluster.TranscriptUtils;
import npTranscript.run.ProcessReadFile.Transcripts;



class ProcessReads{
	
	
	
	int[][] read_count; // [barcode][transcript_index]
	int[][] indices;  // [barcode][transcript_index]
	int[] indices_size;
	int[] barcode_usage;
	
	Map<Integer, Integer>[] forward;
//	Map<Integer, Integer>[] maps; //maps each transcript index to its col position.

	final Transcripts tr;
	final Transcripts barc;
	//final 	List<String> barcodes;

	final Set<String> barc_remainder;
	

	//Barcodes bc;
	PrintWriter leftover;
	final int id_index;
	final int read_id_index, chrom_index, start_index, end_index, breaks_index;
//	final int bc_index;
	 int bc_str_index;
	 final int bc_type_index;
	final int conf_index;
	
	final int ncols; //number of transcripts
	final int nrows; //number of barcodes
	
//	int remainder;
	
	//private final File remainderFile;
	int readsInRemainderFile=0;
	private final OutputStream rem;
	
	public int getRemainder() {
		// TODO Auto-generated method stub
		return readsInRemainderFile;
	}
	

	public void run() {
		try{
		String str = "";
		 for(int i=0; ((str = br.readLine())!=null); i++){
			// if(str.startsWith("readID")) continue;
			process(str);
		}
		//br.close();
		}catch(Exception exc){
			exc.printStackTrace();
		}finally{
			leftover.flush();
		}
		
	}
	BufferedReader br;
	ProcessReads(Transcripts tr, int ncols,int n_barcodes,   InputStream inp, OutputStream rem, int index) throws IOException{
		firstPass = (index==0);
		this.tr = tr;
		this.barc = new Transcripts(false);
		this.barc_remainder = new HashSet<String>();

		this.rem = rem;
		br = new BufferedReader(new InputStreamReader(inp));
		String head = br.readLine();
		this.prefix=index;
	
		this.nrows = n_barcodes;
		 this.ncols = ncols;
	//	remainder = nrows;
		this.read_count = new int[nrows][ncols];
		this.indices = new int[nrows][ncols];
		this.indices_size = new int[nrows];
		this.barcode_usage = new int[nrows];
		Arrays.fill(indices_size, 0);
		this.forward = new Map[nrows];
		/*this.maps = new Map[nrows];*/
		for(int i=0; i<nrows; i++){
			forward[i] = new HashMap<Integer, Integer>();
			Arrays.fill(indices[i], -1);
		}
	
	//	remainderFile.deleteOnExit();
		this.leftover = new PrintWriter(new OutputStreamWriter(rem));
		
		
	
		List<String >header = Arrays.asList( head.split("\t"));
		this.id_index = header.indexOf("id");
		this.read_id_index =  header.indexOf("readID") ;
		this.chrom_index =  header.indexOf("chrom") ;
		this.start_index=  header.indexOf("startPos") ;
		this.end_index =  header.indexOf("endPos") ;
		this.breaks_index =  header.indexOf("breaks_ref") ;

	//	inds = new int[] {read_id_index, chrom_index, start_index, end_index};
		//this.bc_index = header.indexOf("barcode_index");
		this.bc_str_index = header.indexOf("barcode");
		if(bc_str_index<0) bc_str_index = header.indexOf("source");
		this.conf_index = header.indexOf("confidence");
		this.bc_type_index = header.indexOf("barcode_type");
		this.leftover.print(header.get(id_index)+"\t"+header.get(bc_str_index));
		if(bc_type_index>=0)leftover.print("\t"+header.get(bc_type_index)+"\t" +header.get(conf_index));
		this.leftover.println();
		//rem.
	}
	
//	boolean check = true;
	
	
	
	 int max_col=0;
	
	 boolean firstPass = false;
	// final int[] inds;
	public void  process(String str){
		
		String[] line = str.split("\t");
		String transcript = line[id_index];
		if(ProcessReadFile.chrom_only) transcript = transcript.substring(0,transcript.indexOf("/"));
		else{
			if(ProcessReadFile.truncate) transcript = transcript.substring(0,transcript.lastIndexOf("/"));	
			if(!ProcessReadFile.incl_strand) transcript = transcript.replaceAll("/[+\\-]", "");
		}
		Integer trans_ind = this.tr.getAndAdd(transcript,line, read_id_index, chrom_index, start_index, end_index, breaks_index);
				
		if(firstPass){
			this.tr.increment(trans_ind);
		}
		if(conf_index>=0){
			if(!line[conf_index].equals("high")){
				return;
			}
		}
		if(line[bc_str_index].equals("null")){
		//	System.err.println(line[bc_index]+" "+line[this.bc_str_index]); //NA null
			return;// remainder;
		}
		total_reads++;
		String barcode_;
		if(bc_type_index>=0){
			if(line[bc_type_index].endsWith("revC")){
				barcode_ = SequenceUtil.reverseComplement(line[bc_str_index]);
			}else{
				barcode_ = line[bc_str_index];
			}
		}else{
			barcode_ = line[bc_str_index];

		}
		
		Integer  barc_ind = this.barc.get(barcode_);
		
		
		boolean wasNull=false;
		if(barc_ind==null){
			wasNull = true;
			barc_ind = this.barc.getNext(line[bc_str_index],  false);
		}
		if(barc_ind<indices_size.length && this.indices_size[barc_ind]<ncols){
			if(wasNull) this.barc.add(line[bc_str_index], barc_ind); //only add if we using it
			Integer col_index = this.forward[barc_ind].get(trans_ind);
			if(col_index==null){
				col_index = this.indices_size[barc_ind];
				this.forward[barc_ind].put(trans_ind,  col_index);
				this.indices[barc_ind][col_index] = trans_ind;
				this.indices_size[barc_ind]++;
			}
			this.read_count[barc_ind][col_index]++;
			this.barcode_usage[barc_ind]++;
			if(col_index>max_col) {
				max_col=col_index;
			}
		}else{
			/*if(barc_ind<indices_size.length && indices_size[barc_ind]==ncols){
				remainder--;
				indices_size[barc_ind]++;
			}*/
			this.leftover.println(line[id_index]+"\t"+line[bc_str_index]+(conf_index<0 ? "": "\t"+line[bc_type_index] +"\t"+line[conf_index]));
			this.readsInRemainderFile++;
			this.barc_remainder.add(line[bc_str_index]);
		}
		//return remainder;
//		this.cluster_ids.in
	}
	final int prefix;
	public int total_reads=0;
	
	public int remainingBarcodes(){
		return barc_remainder.size();
	}
	
	static class IntInt implements Comparable {
		int index; int count;
		
		@Override
		public int compareTo(Object o) {
			return Integer.compare(((IntInt)o).count,count);
					
		}
	}
	static class IntIntC{
		final IntInt[] obj;
		IntIntC(int len){
			obj = new IntInt[len];
			for(int i=0; i<len; i++){
				obj[i] = new IntInt();
			}
		}
		int[] counts;
		int[] indices;
		int len;
		void fill(int[] counts, int[] indices){
			this.counts = counts;
			this.indices = indices;
			this.len = counts.length;
			for(int i=0; i<len; i++){
				obj[i].count = counts[i];
				obj[i].index = indices[i];
			}
		}
		
		public void fill(int[] counts) {
			this.len = counts.length;
			this.counts = counts;
			indices = new int[len];
			for(int i=0; i<len; i++){
				indices[i] = i;
				obj[i].count = counts[i];
				obj[i].index = i;
			}
		}
		
		void sort(){
			Arrays.sort(obj);
			for(int i=0;i<obj.length; i++){
				counts[i] = obj[i].count;
				indices[i] = obj[i].index;
			}
		}
		
	}
	
	static int[][] resizeColsAndRows(int[][] read_count, int len1, int ncol1){
		int[][]read_count1 = new int[len1][ncol1];
		for(int i=0; i<len1;i++){
			System.arraycopy(read_count[i], 0, read_count1[i], 0, ncol1);
		}
		return read_count1;
	}
	static int[][] resizeRows(int[][] read_count, int len1){
		int[][]read_count1 = new int[len1][];
		System.arraycopy(read_count, 0, read_count1, 0, len1);
		return read_count1;
	}
	static int[] resizeRows(int[] read_count, int len1){
		int[]read_count1 = new int[len1];
		System.arraycopy(read_count, 0, read_count1, 0, len1);
		return read_count1;
	}
	private static  int[][] reorderRows(int[][] read_count, int[] indices) {
		int len2 = indices.length;
		int[][] counts2 = new int[len2][];
		for(int i=0; i<len2; i++){
			counts2[i] = read_count[indices[i]];
		}
		return counts2;
	}
	private static  String[] reorderRows(String[] read_count, int[] indices) {
		int len2 = indices.length;
		String[] counts2 = new String[len2];
		for(int i=0; i<len2; i++){
			counts2[i] = read_count[indices[i]];
		}
		return counts2;
	}
	private static int[] reorderRows(int[] read_count, int[] indices) {
		int len2 = indices.length;
		int[] counts2 = new int[len2];
		for(int i=0; i<len2; i++){
			counts2[i] = read_count[indices[i]];
		}
		return counts2;
	}
	
	public void finalise(IHDF5Writer altT, boolean reorder_col, boolean reorder_row) throws IOException{
		//leftover.close();
		int barc_len = this.indices_size.length;
		int len1 = barc_len;
		while(this.indices_size[len1-1]==0){
			len1 = len1-1;
		}
		String[] barc_array = this.barc.getArray(0,len1);
		if(max_col<this.ncols-1){
			//make col_size smaller to fit data also fit to max cols
			int ncol1 = max_col+1;
			System.err.println("less transcript cols than expected "+ncol1+" vs "+this.ncols);
			read_count = resizeColsAndRows(read_count, len1, ncol1);
			indices = resizeColsAndRows(indices, len1, ncol1);
			if(len1<nrows){
				barcode_usage =resizeRows(barcode_usage,len1);
			}
		//	System.err.println("h");
		}else if(len1 < this.nrows){
			System.err.println("less rows than expected "+len1+" vs "+this.nrows);
			read_count = resizeRows(read_count, len1);
			indices = resizeRows(indices,len1);
			barcode_usage =resizeRows(barcode_usage,len1);
		}
		
		
		if(reorder_row){
			IntIntC  obj = new IntIntC(barcode_usage.length);
			obj.fill(barcode_usage);
			obj.sort(); // this sounds the barcode_usage, so we dont resort it
			read_count = reorderRows(read_count, obj.indices);
			indices = reorderRows(indices, obj.indices);
			barc_array = reorderRows(barc_array, obj.indices);
			//barcode_usage = reorderRows(barcode_usage, obj.indices); /// dont resort
		}
		if(reorder_col){
			IntIntC  obj = new IntIntC(read_count[0].length);
			for(int i=0; i<read_count.length; i++){
				obj.fill(read_count[i], indices[i]);
				obj.sort();
			}
		}
		
		altT.writeIntMatrix(prefix+"/counts", read_count);
		altT.writeIntMatrix(prefix+"/indices", indices);
		altT.writeIntArray(prefix+"/barcode_usage", barcode_usage);
		altT.writeStringArray(prefix+"/barcodes", barc_array);
		//altT.close();
	}


	

	
	
}