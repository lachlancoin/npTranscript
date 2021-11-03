package npTranscript.run;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FilenameFilter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.OutputStream;
import java.io.OutputStreamWriter;
import java.io.PrintWriter;
import java.io.SequenceInputStream;
import java.util.ArrayList;
import java.util.Enumeration;
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

	//public static String barcode_file=null;
	//public static String leftOverFile = null;
//	private static final Logger LOG = LoggerFactory.getLogger(HTSErrorAnalysisCmd.class);

	public ProcessReadFile() {
		super();
		Deployable annotation = getClass().getAnnotation(Deployable.class);
		setUsage(annotation.scriptName() + " [options]");
		setDesc(annotation.scriptDesc());
		addBoolean("overwrite", true, "whether to delete and start from scratch");
		addBoolean("include_strand", true, "whether to remove strand from barcode");
		addBoolean("truncate", true, "whether to only include locus and not isoform information ");


		addBoolean("reorder", true, "whether to re-order barcodes and transcripts");
	//	addBoolean("saveReadID", true, "whether to save the readID of first read in each cluster");

		addString("barcode_file", null, "Name of barcode file", false);
		addString("inputFileDir", null, "Name of input file", true);
		addString("outputFile", "out.h5", "Name of output file", false);
		addString("refFileDir", null, "reference input", true);
		
		addString("suffix", ".reads.txt.gz", "suffix of input files", false);
		
		this.addDouble("maxcells", 4e6, "Max number of cells per block",false); 
		this.addInt("num_barcodes", 4000, "Max number of barcodes",false);

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
	 static boolean truncate;
	 static boolean incl_strand;
	 static boolean overwrite;
	 static double maxcells;
	 //static int num_barcodes;
	 static boolean reorder;
	 static String suffix;
	 
	 static  FilenameFilter filter = new FilenameFilter(){

			@Override
			public boolean accept(File dir, String name) {
				return name.endsWith(suffix);
			}
			 
	 };
	public static void main(String[] args1) throws IOException, InterruptedException {
		tme0 = System.currentTimeMillis();
		CommandLine cmdLine = new ProcessReadFile();
		String[] args = cmdLine.stdParseLine(args1);
		String[] inpDirs =  cmdLine.getStringVal("inputFileDir").split(":");
		
		File refDir = new File( cmdLine.getStringVal("refFileDir"));

		int num_barcodes=cmdLine.getIntVal("num_barcodes");
		 reorder = cmdLine.getBooleanVal("reorder");
		 maxcells = cmdLine.getDoubleVal("maxcells");
		truncate = cmdLine.getBooleanVal("truncate");
		incl_strand = cmdLine.getBooleanVal("include_strand");
		overwrite=cmdLine.getBooleanVal("overwrite");
		suffix = cmdLine.getStringVal("suffix");
		
		String  outFile = cmdLine.getStringVal("outputFile");		
		System.err.println("doing reference");
		File transcripts_file =  run(refDir, outFile,  null,1); //only one barcode for ref
		 System.err.println("doing main class");
		 for(int i=0 ; i<inpDirs.length; i++){
			 System.err.println("running "+i);
			 run(new File(inpDirs[i]),  outFile, transcripts_file, num_barcodes);
		 }
	}
	
	static InputStream getInputStream(File dir, final String[] file){
		if(file==null) return System.in;
		Enumeration<InputStream> en = new Enumeration<InputStream>(){
			int i=0;
			@Override
			public boolean hasMoreElements() {
				return i < file.length;
			}

			@Override
			public InputStream nextElement() {
				InputStream inp =  null;
				try{
				inp = new FileInputStream(new File(dir, file[i]));
				if(file[i].endsWith(".gz")) inp = new GZIPInputStream(inp);
				}catch(IOException exc){
					exc.printStackTrace();
				}
				// TODO Auto-generated method stub
				i++;
				
				return inp;
			}
			
		};
		return new SequenceInputStream(en);
	}
	
	public static File  run(File inpDir, String outFileName, File transcripts_file, int num_barcodes1) throws IOException{
		File outFile = new File(inpDir, outFileName);
		if(overwrite) outFile.delete();
		String[] input_file = inpDir.list(filter);
		IHDF5Writer altT = HDF5Factory.open(outFile);
		//altT.writeStringArray("/barcodes", barcodes.toArray(new String[0]));//.subList(0, len1).toArray(new String[0]));

		InputStream inp = getInputStream(inpDir, input_file);
		int remainder = Integer.MAX_VALUE;
		
		File transcripts_out_file =new File( outFile.getAbsolutePath()+".transcripts.txt.gz");
		
		File transcripts_out_file1 = new File(outFile.getAbsolutePath()+".transcripts.mod.txt.gz");

		PrintWriter transcripts_pw = new PrintWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(transcripts_out_file))));
		boolean addCounts = true;
		
		Transcripts tr = new Transcripts(transcripts_pw, transcripts_file, addCounts);
		for(int index=0; remainder>0;index++){
			int  ncols = (int) Math.floor(maxcells/(double)num_barcodes1);
			
			if(remainder < ncols) ncols = remainder; // worst case if each read is different cluster
			double cells = ncols *num_barcodes1;
			System.err.println("ncols "+ncols+" by  "+num_barcodes1+" "+String.format("%5.3g", cells));
			File remainderFile = new File("tmp.leftover."+System.currentTimeMillis()+"."+Math.random());
			remainderFile.deleteOnExit();
			OutputStream rem = new FileOutputStream(remainderFile);
			ProcessReads pr = 
					new ProcessReads(tr, ncols,num_barcodes1, 	inp,rem,  index);
			pr.run();
			pr.finalise(altT, reorder, reorder);
			if(index==0){
				transcripts_pw.close();
				if(addCounts){
					tr.appendCounts(transcripts_out_file, transcripts_out_file1, 100);
					transcripts_out_file.delete();
				}
			}
			inp.close();
			rem.close();
			remainder=pr.getRemainder();
			if(remainder>0){
			inp =new FileInputStream(remainderFile);
			}else{
				inp=null;
			}
			num_barcodes1 = Math.min(num_barcodes1, pr.remainingBarcodes());
			System.err.println(num_barcodes1+" barcodes_remaining");
			 System.err.println("reads remaining "+remainder+ " of total "+tr.cluster_ids.size());
			
		}
	//	altT.writeStringArray("/transcripts", tr.getArray()); //need to consider the offset
		
		altT.close();
		return addCounts ? transcripts_out_file1: transcripts_out_file;
	}
	
	// keeps a uniform list of transcripts
	public static class Transcripts{
		//final boolean saveReadID;
		public Transcripts(boolean addCounts){
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

		public Transcripts(PrintWriter pw, File inputF, boolean addCounts) {
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
		private List<Integer> counts = null;//new ArrayList<Integer>();
		
		private Map<String, Integer> cluster_id_map = new HashMap<String, Integer>();
		public Integer get(String transcript){
			return this.cluster_id_map.get(transcript);
		}
		
		 boolean addcount;
		//assumes is absent
		public Integer getNext(String transcript,  boolean add){
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
		public Integer getAndAdd(String transcript, String[] line,
				int read_id_index, int chrom_index, int start_index, int end_index
				) {

			
			Integer trans_ind = this.cluster_id_map.get(transcript);
			if(trans_ind==null){
				trans_ind = getNext(transcript,  true);
				if(pw!=null){
					pw.println(transcript+"\t"+
							line[read_id_index]+"\t"+
							line[chrom_index]+"\t"+ line[start_index]+"\t"+line[end_index]
							);
				}
			}
			
			return trans_ind;
		}
		public void add(String transcript, Integer trans_ind){
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
				res[i] = this.cluster_ids.get(i+offset);
			}
			return res;
		}
	}
	
	
	
	
	
}