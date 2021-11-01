package npTranscript.run;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.OutputStream;
import java.io.OutputStreamWriter;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;

import ch.systemsx.cisd.hdf5.HDF5Factory;
import ch.systemsx.cisd.hdf5.IHDF5Writer;
import japsa.util.CommandLine;
import japsa.util.deploy.Deployable;

@Deployable(scriptName = "npTranscript.run", scriptDesc = "Analysis of coronavirus sequence data")
public class ProcessReadFile extends CommandLine {

	public static String barcode_file=null;
	public static String leftOverFile = null;
public static String prefix="";
//	private static final Logger LOG = LoggerFactory.getLogger(HTSErrorAnalysisCmd.class);

	public ProcessReadFile() {
		super();
		Deployable annotation = getClass().getAnnotation(Deployable.class);
		setUsage(annotation.scriptName() + " [options]");
		setDesc(annotation.scriptDesc());
		addBoolean("overwrite", false, "whether to delete and start from scratch");
		addString("barcode_file", null, "Name of barcode file", false);
		addString("prefix", "0", "Name of prefix", false);
		addString("inputFile", null, "Name of input file", false);
		addString("leftOverFile",null,"Name of file for leftover, or stdout");
		this.addInt("ncols", 100, "Max number of transcripts",false);

	}
 static 	long tme0;
	public static void main(String[] args0) throws IOException, InterruptedException {
		tme0 = System.currentTimeMillis();
		boolean stdin = false;
		String[] args1 = args0;
		if(args0[args0.length-1].equals("-")){
			stdin = true;
			args1 = new String[args0.length-1];
			System.arraycopy(args0, 0, args1, 0, args1.length);
		}
		CommandLine cmdLine = new ProcessReadFile();
		String[] args = cmdLine.stdParseLine(args1);
		String barcode_files=cmdLine.getStringVal("barcode_file");
		String input_file = cmdLine.getStringVal("inputFile");
		prefix = cmdLine.getStringVal("prefix");
		leftOverFile = cmdLine.getStringVal("leftOverFile");
		List<String> barcodes = new ArrayList<String>();
		{
			InputStream is = new FileInputStream(barcode_files);
			if(barcode_files.endsWith(".gz")) is = new GZIPInputStream(is);
			BufferedReader br = new BufferedReader(new InputStreamReader(is));
			String st = "";
			while((st = br.readLine())!=null){
				barcodes.add(st.split("-")[0]);
			}
			br.close();
		}
			System.err.println("using barcodes"+Arrays.asList(barcode_files));
		Integer ncols = cmdLine.getIntVal("ncols");
	
		

		InputStream inp =new FileInputStream(input_file);
		if(input_file.endsWith(".gz")) inp = new GZIPInputStream(inp);
		OutputStream leftover = leftOverFile==null ? System.out : new FileOutputStream(leftOverFile);
		if(leftOverFile!=null && leftOverFile.endsWith(".gz")) leftover = new GZIPOutputStream(leftover);
		boolean truncate = false;
		
		BufferedReader br = new BufferedReader(new InputStreamReader(inp));
		String head = br.readLine();
		ProcessReads pr = 
				new ProcessReads( ncols, barcodes, 	leftover, truncate, head);
		String str = "";
		int max_reads =Integer.MAX_VALUE;
		try{
		inner: for(int i=0; i<max_reads &&((str = br.readLine())!=null); i++){
		
			int remainder = pr.process(str);
			if(remainder<=0){
				System.err.println("is full");
				break inner;
			}
			
		}
		}catch(Exception exc){
			exc.printStackTrace();
		}
	
		boolean overwrite=cmdLine.getBooleanVal("overwrite");
		File outFile = new File("out.h5");

		if(overwrite) outFile.delete();
		IHDF5Writer altT = HDF5Factory.open(outFile);
		altT.writeStringArray("/barcodes", barcodes.toArray(new String[0]));//.subList(0, len1).toArray(new String[0]));
		pr.finalise(altT, prefix);
		altT.close();
		br.close();
		leftover.close();
	}
	

	
static class ProcessReads{
	int[][] read_count; // [barcode][transcript_index]
	int[][] indices;  // [barcode][transcript_index]
	int[] indices_size;
	
	Map<Integer, Integer>[] forward;
//	Map<Integer, Integer>[] maps; //maps each transcript index to its col position.

	List<String>cluster_ids = new ArrayList<String>();
	Map<String, Integer> cluster_id_map = new HashMap<String, Integer>();
	
	List<Integer>barcode_ids = new ArrayList<Integer>();
	Map<Integer, Integer> barcode_id_map = new HashMap<Integer, Integer>();
	

	//Barcodes bc;
	PrintWriter leftover;
	final int id_index;
	final int bc_index;
	final int bc_str_index;
	final int conf_index;
	
	final int ncols; //number of transcripts
	final int nrows; //number of barcodes
	
	int remainder;
final 	List<String> barcodes;
	
	ProcessReads(int ncols,List<String> barcodes, OutputStream leftover, boolean truncate, String head) throws IOException{
		this.truncate  =truncate;
		this.barcodes = barcodes;
		this.nrows = barcodes.size();
		 this.ncols = ncols;
		remainder = nrows;
		this.read_count = new int[nrows][ncols];
		this.indices = new int[nrows][ncols];
		this.indices_size = new int[nrows];
		Arrays.fill(indices_size, 0);
		this.forward = new Map[nrows];
		/*this.maps = new Map[nrows];*/
		for(int i=0; i<nrows; i++){
			forward[i] = new HashMap<Integer, Integer>();
			Arrays.fill(indices[i], -1);
		}
		this.leftover = new PrintWriter(new OutputStreamWriter(leftover));
		
		
		this.leftover.println(head);
		List<String >header = Arrays.asList( head.split("\t"));
		this.id_index = header.indexOf("id");
		this.bc_index = header.indexOf("barcode_index");
		this.bc_str_index = header.indexOf("barcode");
		this.conf_index = header.indexOf("confidence");
	}
	
//	boolean check = true;
	
	
	final boolean truncate;
	
	public int  process(String str){
		String[] line = str.split("\t");
		String transcript = line[id_index];
		if(truncate) transcript = transcript.substring(0,transcript.lastIndexOf("/"));	
		Integer trans_ind = this.cluster_id_map.get(transcript);
		if(trans_ind==null){
			cluster_id_map.put(transcript,  trans_ind=this.cluster_ids.size());
			this.cluster_ids.add(transcript);
		}
		if(line[bc_index].equals("NA")){
			return remainder;
		}
		Integer  barcode = Integer.parseInt(line[bc_index]);
		Integer barc_ind = this.barcode_id_map.get(barcode);
		if(barc_ind==null){
			barcode_id_map.put(barcode, barc_ind = this.barcode_ids.size());
			this.barcode_ids.add(barcode);
		}
		
		if(this.indices_size[barc_ind]<ncols){
		
			Integer col_index = this.forward[barc_ind].get(trans_ind);
			if(col_index==null){
				col_index = this.indices_size[barc_ind];
				this.forward[barc_ind].put(trans_ind,  col_index);
				this.indices[barc_ind][col_index] = trans_ind;
				this.indices_size[barc_ind]++;
			}
			this.read_count[barc_ind][col_index]++;
		}else{
			if(indices_size[barc_ind]==ncols){
				remainder--;
				indices_size[barc_ind]++;
			}
		//	System.err.println(barc_ind+" is leftover");
			this.leftover.println(str);
			
		}
		return remainder;
//		this.cluster_ids.in
	}
	
	public void finalise(IHDF5Writer altT, String prefix) throws IOException{
		leftover.close();
		int barc_len = this.barcode_ids.size();

		int len1 = barc_len;
		while(this.indices_size[len1-1]==0){
			len1 = len1-1;
		}
		

		
		
		int[] barcode_ids_ = new int[len1];
		for(int i=0; i<len1; i++){
			barcode_ids_[i] = this.barcode_ids.get(i);
		}
		altT.writeIntArray(prefix+"/barcode_indices", barcode_ids_);
		altT.writeStringArray(prefix+"/transcripts", this.cluster_ids.toArray(new String[0]));
		int[][] read_count1, indices1;
		if(len1 < this.nrows){
			read_count1 = new int[len1][];
			indices1 = new int[len1][];
			System.arraycopy(read_count, 0, read_count1, 0, len1);
			System.arraycopy(indices, 0, indices1, 0, len1);
		}else{
			read_count1 = this.read_count;
			indices1 = this.indices;
		}
		altT.writeIntMatrix(prefix+"/counts", read_count1);
		altT.writeIntMatrix(prefix+"/indices", indices1);
		//altT.close();
	}
	
}
	
	
	
}