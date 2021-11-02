package npTranscript.run;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.OutputStream;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.zip.GZIPInputStream;

import ch.systemsx.cisd.hdf5.HDF5Factory;
import ch.systemsx.cisd.hdf5.IHDF5Writer;
import japsa.util.CommandLine;
import japsa.util.deploy.Deployable;

@Deployable(scriptName = "npTranscript.run", scriptDesc = "Analysis of coronavirus sequence data")
public class ProcessReadFile extends CommandLine {

	//public static String barcode_file=null;
	//public static String leftOverFile = null;
//	private static final Logger LOG = LoggerFactory.getLogger(HTSErrorAnalysisCmd.class);

	public ProcessReadFile() {
		super();
		Deployable annotation = getClass().getAnnotation(Deployable.class);
		setUsage(annotation.scriptName() + " [options]");
		setDesc(annotation.scriptDesc());
		addBoolean("overwrite", false, "whether to delete and start from scratch");
		addBoolean("reorder", true, "whether to re-order barcodes and transcripts");
	//	addBoolean("saveReadID", true, "whether to save the readID of first read in each cluster");

		addString("barcode_file", null, "Name of barcode file", false);
		addString("inputFile", null, "Name of input file", false);
		addString("outputFile", "out.h5", "Name of output file", false);
		addString("ref_transcripts", null, "Annotated list of transcripts to use", false);
		this.addDouble("maxcells", 1e6, "Max number of cells per block",false); 
		this.addInt("num_barcodes", 2000, "Max number of barcodes",false);

	}
 static 	long tme0;
 
 
 
 public static int countBarcode(String barcode_files) throws IOException{
	 int num_barcodes=0;
			InputStream is = new FileInputStream(barcode_files);
			if(barcode_files.endsWith(".gz")) is = new GZIPInputStream(is);
			BufferedReader br = new BufferedReader(new InputStreamReader(is));
			String st = "";
			while((st = br.readLine())!=null){
				num_barcodes++;
			}
			br.close();
			return num_barcodes;
 }
 
	public static void main(String[] args1) throws IOException, InterruptedException {
		tme0 = System.currentTimeMillis();
		CommandLine cmdLine = new ProcessReadFile();
		String[] args = cmdLine.stdParseLine(args1);
		//String barcode_files=cmdLine.getStringVal("barcode_file");
		String input_file = cmdLine.getStringVal("inputFile");
		//prefix = cmdLine.getStringVal("prefix");
		//leftOverFile = cmdLine.getStringVal("leftOverFile");
		Integer num_barcodes=cmdLine.getIntVal("num_barcodes");
		//if(barcode_files!=null) num_barcodes = countBarcode(barcode_files);
	//	List<String> barcodes = new ArrayList<String>();
		boolean reorder = cmdLine.getBooleanVal("reorder");
	//	boolean saveReadID = cmdLine.getBooleanVal("saveReadID");
	//		System.err.println("using barcodes"+Arrays.asList(barcode_files));
		Double maxcells = cmdLine.getDoubleVal("maxcells");
		
	//	int max_reads =Integer.MAX_VALUE;
		
		File outFile = new File(cmdLine.getStringVal("outputFile"));
		boolean truncate = false;
		boolean overwrite=cmdLine.getBooleanVal("overwrite");
		if(overwrite) outFile.delete();
		IHDF5Writer altT = HDF5Factory.open(outFile);
		//altT.writeStringArray("/barcodes", barcodes.toArray(new String[0]));//.subList(0, len1).toArray(new String[0]));

		InputStream inp =input_file==null ? System.in : new FileInputStream(input_file);
		if(input_file!=null && input_file.endsWith(".gz")) inp = new GZIPInputStream(inp);
		int remainder = Integer.MAX_VALUE;
		String transcripts_file =cmdLine.getStringVal("ref_transcripts");
		PrintWriter transcripts_pw = new PrintWriter(new FileWriter(outFile+".transcripts.txt"));

		Transcripts tr = new Transcripts(transcripts_pw, transcripts_file);
		for(int index=0; remainder>0;index++){
			int  ncols = (int) Math.floor(maxcells/num_barcodes.doubleValue());
			System.err.println("reads remaining "+remainder);
			if(remainder < ncols) ncols = remainder; // worst case if each read is different cluster
			double cells = ncols *num_barcodes;
			System.err.println("ncols "+ncols+" by  "+num_barcodes+" "+String.format("%5.3g", cells));
			File remainderFile = new File("tmp.leftover."+System.currentTimeMillis()+"."+Math.random());
			remainderFile.deleteOnExit();
			OutputStream rem = new FileOutputStream(remainderFile);
			ProcessReads pr = 
					new ProcessReads(tr, ncols,num_barcodes, 	truncate, inp,rem,  index);
			pr.run();
			pr.finalise(altT, reorder, reorder);
			inp.close();
			rem.close();
			remainder=pr.getRemainder();
			if(remainder>0){
			inp =new FileInputStream(remainderFile);
			}else{
				inp=null;
			}
			num_barcodes = Math.min(num_barcodes, pr.remainingBarcodes());
			System.err.println(num_barcodes+" barcodes_remaining");
			
			
		}
		transcripts_pw.close();
	//	altT.writeStringArray("/transcripts", tr.getArray()); //need to consider the offset
		
		altT.close();
		
	}
	
	// keeps a uniform list of transcripts
	public static class Transcripts{
		//final boolean saveReadID;
		public Transcripts(){
			
		}
		
		public Transcripts(PrintWriter pw, String input) {
			//this.saveReadID = saveReadID;
			this.pw = pw;
			if(input!=null){
				try{
					File inputF = new File(input);
					if(inputF.exists()){
				BufferedReader br = new BufferedReader(new FileReader(inputF));
				String str = "";
				for(int i=0; (str=br.readLine())!=null; i++){
					String[] st = str.split("\t");
					cluster_id_map.put(st[0], i);
					this.cluster_ids.add(st[0]);
					pw.println(str);
				}
				br.close();
					}
				}catch(IOException exc){
					exc.printStackTrace();
				}
			}
		}
		
		PrintWriter pw = null;
		private List<String>cluster_ids = new ArrayList<String>();
		private Map<String, Integer> cluster_id_map = new HashMap<String, Integer>();
		public Integer get(String transcript){
			return this.cluster_id_map.get(transcript);
		}
		//assumes is absent
		public Integer getNext(String transcript, String readID, boolean add){
			Integer trans_ind = cluster_ids.size();
			if(add){
				cluster_id_map.put(transcript,  trans_ind);
				this.cluster_ids.add(transcript);
			
			}
			return trans_ind;
		}
		public Integer getAndAdd(String transcript, String readID) {
			Integer trans_ind = this.cluster_id_map.get(transcript);
			if(trans_ind==null){
				trans_ind = getNext(transcript, readID, true);
				if(pw!=null){
					pw.print(transcript);
					pw.print("\t"+ readID);
					pw.println();
				}
			}
			return trans_ind;
		}
		public void add(String transcript, Integer trans_ind){
			cluster_id_map.put(transcript,  trans_ind);
			this.cluster_ids.add(transcript);
		}
		public String[] getArray() {
			return cluster_ids.toArray(new String[0]);
		}
		public String[] getArray(int offset, int len1) {
			String[] res = new String[len1];
			for(int i=0; i<len1; i++){
				res[i] = this.cluster_ids.get(i+offset);
			}
			return res;
		}
	}
	
	
	
	
	
}