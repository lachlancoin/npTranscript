/*****************************************************************************
 * Copyright (c)Lachlan Coin University of Melbourne, All rights reserved.         *
 *                                                                           *
 * Redistribution and use in source and binary forms, with or without        *
 * modification, are permitted provided that the following conditions        *
 * are met:                                                                  * 
 *                                                                           *
 * 1. Redistributions of source code must retain the above copyright notice, *
 *    this list of conditions and the following disclaimer.                  *
 * 2. Redistributions in binary form must reproduce the above copyright      *
 *    notice, this list of conditions and the following disclaimer in the    *
 *    documentation and/or other materials provided with the distribution.   *
 * 3. Neither the names of the institutions nor the names of the contributors*
 *    may be used to endorse or promote products derived from this software  *
 *    without specific prior written permission.                             *
 *                                                                           *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS   *
 * IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, *
 * THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR    *
 * PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR         *
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,     *
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,       *
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR        *
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF    *
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING      *
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS        *
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.              *
 ****************************************************************************/

/*                                                    
 * 01/03/2020 - Lachlan Coin                                       
 ****************************************************************************/

package npTranscript.run;


import java.io.BufferedReader;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintStream;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Date;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;
import java.util.concurrent.Executors;
import java.util.regex.Pattern;
import java.util.stream.Collectors;
import java.util.zip.GZIPOutputStream;

import com.google.gson.Gson;

import htsjdk.samtools.SAMProgramRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import japsa.seq.Alphabet;
import japsa.seq.Sequence;
import japsa.seq.SequenceReader;
import japsa.tools.seq.SequenceUtils;
import japsa.util.CommandLine;
import japsa.util.deploy.Deployable;
import npTranscript.cluster.Annotation;
import npTranscript.cluster.CigarCluster;
import npTranscript.cluster.CigarHash;
import npTranscript.cluster.CigarHash2;
import npTranscript.cluster.EmptyAnnotation;
import npTranscript.cluster.GFFAnnotation;
import npTranscript.cluster.IdentityProfile1;
import npTranscript.cluster.IdentityProfileHolder;
import npTranscript.cluster.Outputs;
import npTranscript.cluster.TranscriptUtils;

/**
 * @author Lachlan Coin
 *
 */

@Deployable(scriptName = "npTranscript.run", scriptDesc = "Analysis of coronavirus sequence data")
public class ViralTranscriptAnalysisCmd2 extends CommandLine {

	
	static class SAMRecord1  implements Comparable{
		SAMRecord1(SAMRecord sr){
			this.sr = sr;
		}
		SAMRecord sr;
		@Override
		public int compareTo(Object o) {
			SAMRecord1 s1 = (SAMRecord1)o;
		    int st = this.sr.getAlignmentStart();
		    int st1 = s1.sr.getAlignmentStart();
		    if(st <st1) return -1;
		    if(st>st1) return 1;
		    else return 0;
			   
		}
		
	}
	public static boolean exclude_indeterminate_strand=true;
	public static boolean exclude_polyT_strand=false;
	public static boolean exclude_reads_without_barcode=false;
	public static boolean verbose=false;
	public static boolean illumina =false;
	public static boolean reverseForward = true;
	public static int max_umi = 100; // including umi + barcode.  This should be bit bigger for getting polyT length , i.e. if reverseForward=false
	public static int min_umi = 15;// including umi + barcode
	
public static int readsToSkip=0;
	

//	private static final Logger LOG = LoggerFactory.getLogger(HTSErrorAnalysisCmd.class);

	public ViralTranscriptAnalysisCmd2() {
		super();
		Deployable annotation = getClass().getAnnotation(Deployable.class);
		setUsage(annotation.scriptName() + " [options]");
		setDesc(annotation.scriptDesc());
		
		addString("readsOutputFile", null, "Name of file for output", false);
		addString("flames_barcodes", null, "Extracted flames barcdes", false);
		addString("optsFile", null, "Name of file with extra options", false);
		addBoolean("verbose", false, "verbose output", false);
		addString("types_to_include",null, "For coronavirus only, which types to include, e.g. 5_3:5_no3", false);
		addString("optsType", null, "Which column in opts file", false);
		addInt("max_umi",150,"Maximum size of UMI + barcode and possibly polyT minus bc_len_AT if reverse forward is true");
		addInt("min_umi",15,"Minimum size of UMI + barcode");
		addInt("tolerance_barcode",1,"Maximum edit dist to barcode");
		addBoolean("reverseForward", false, "For finding polyA tail whether to reverse the polyT before edlib step.  If false will include some polyT ");
		addInt("barcode_extent", 200, "search for barcode in first Xbp");
		addInt("barcode_ignore", 0, "search for barcode in first Ybp");
		addString("barcode_list",null, "list for decoding barcodes (used if more than one bamfile when streaming from fastq");
		addString("barcode_file",null, "barcode file");
		addString("annotation_mode","none", "annotation mode");
		

		addString("api_url",null, "URL of database");
		addInt("report",5000,"number of reads to aggregate in single line of json output");	
		
/*		addInt("bc_len_AT",6,"Length of polyA or polyT tract to look for at ends");
		addInt("edit_thresh_AT",1,"Max edit dist (inclusive for AT tracts");
		addInt("max_dist_to_end_AT",20,"Max distance from end of reads to find AT");
		addInt("dist_to_ignore_AT",0,"Distance to ignore for polyA at ends of read");
*/

		addBoolean("exclude_indeterminate_strand", false,"whether to exclude if we cnat figure out strand for cDNA");
		addBoolean("exclude_polyT_strand", false,"whether to exclude RNA if it has a polyT strand");
		addBoolean("exclude_reads_without_barcode", false,"whether to exclude read if it does not have barcode");
addBoolean("illumina", false, "use illumina libary");
		
		addString("inputFile", "-", "Name of input file", false);
		addString("sampleID",null, "sampleID");
		//addString("fastqFile", null, "Name of bam file", false);
		addString("reference", null, "Name of reference genome", false);
		addString("annotation_file", null, "ORF annotation file or GFF file", false);
		addBoolean("useExons", true, "wehether to use exons");
	//	addBoolean("baseBreakPointsOnFirst", false, "whether the break points should be preferentially based on the first bam");
		addBoolean("sorted", false, "whether bamfile sorted in terms of position.");
		addBoolean("sequential", true, "whether to iterate through one bam at a time (in sequence) or in parallel.");
	
		addBoolean("enforceStrand", false, "whether to enforce strandedness");
		addString("readList", "", "List of reads", false);
		addString("join_out",null, "Output for joining sequences");
		addString("overlap_out",null, "Output for overlap sequences");

		addString("chroms_to_include", "all", "Restrict to these chroms, colon delimited", false);
		addString("chroms_to_ignore", "none", "Ignore these chroms", false);
	//	addString("bedChr", null, "Use this for the chrom in bed chr, e.g. NC_045512v2, false");
        addBoolean("singleGFF", false, "whether to have single output gff, rather than different for different input. ");
	//	addString("resdir", "results"+System.currentTimeMillis(), "results directory");
		addString("RNA", "name", "If is direct RNA.  Can be tab delimmited boolean, e.g. true:true  or if set to name it will look for RNA string in name");
		addBoolean("trainStrand", false,"whether to produce training data for strand correction");
		addBoolean("recordDepthByPosition", false, "whether to store position specific depth (high memory");
	//	addString("annotToRemoveFromGFF",null, "annotation to remove from GFF , BED and ref files");
		addInt("maxReads", Integer.MAX_VALUE, "max no of reads");
		addInt("readsToSkip",0, "no_reads_to_skip");
		//addBoolean("")
		addDouble("probInclude", 1.0, "probability of including each read");
		addInt("minClusterEntries",10,"threshold for consensus");
	//	addBoolean("tryComplementOnExtra", false, "look for negative strand matches on left over seqs");
		addBoolean("reAlignExtra", false, "whether to try realigning the extra sequence");
	//	addBoolean("combineOutput", false, "whether to combine output from different chroms");
		addString("pattern", null, "Pattern of read name, used for filtering");
		addString("span", "protein_coding", "Filtering span.  Use all not to filter.");
		addInt("qual", 0, "Minimum quality required");
		
		addString("round", "1", "Bin size for numerical hashing");
	//	addInt("bin1",100, "Bin size for hashing of the end");
		addDouble("overlap_thresh",0.33,"max percentage overlap between supp alignments considered");
		addInt("overlap_max",20,"max bp overlap on read between supp alignments considered");
		addInt("gap_max",20,"min bp gap on read between supp alignments considered");

		addString("numExonsMSA", "none", "number of exons For  MSA");
		addInt("max_seqs_per_cluster",100000,"max sequences in MSA files");
		addInt("breakThresh", 1000, "Thresh for break points to match clusters.  If bigger than genome size then no break points");
		addBoolean("includeStart", true, "Whether to include start position in the cluster hash");
		addString("chromsToRemap", "", "names of chromosomes (entries in reference) which should be remapped (colon delimited) .  Will produce fastq output files for each chrom with mapping reads");
		addInt("coverageDepthThresh", 100, "Threshhold for writing base level depth information to h5 file");
		addString("isoformDepthThresh", "10", "Threshhold for printing out all isoforms");
		addDouble("msaDepthThresh", 10, "Threshhold for running MSA per subcluster");
		addDouble("leftOverQualThresh", 20, "Quality thresh for leftover seqs");
		addString("fail_thresh", "0:0", "Pass threshold (first all reads and second for reads which do not contain any existing annotation");
		addString("doMSA", "false" , "Options: 5_3 or all or span=0 ");
		addString("library",".","library for h5 files");
		addString("msa_source", null , "indicates how to group msa, e.g 0,1,2;3,4,5 ");
		addString("coords_msa","coords.txt", "which coords to keep");
		//addString("aligner", "kalign" , "Options: kalign3, poa, spoa, abpoa");
		addInt("startThresh", 100, "Threshold for having 5'");
		addInt("endThresh", 100, "Threshold for having 3'");
		addInt("maxThreads", 8, "max threads for processing");
		addBoolean("calcBreaks", false, "whether to calculate breakpoints");
		addInt("extra_threshold", 200, "Threshold saving umatched 3'or 5'parts of reads");
		addInt("min_first_last_exon",0, "minimum first or last exon size");
		//addDouble("overlapThresh", 0.95, "Threshold for overlapping clusters");
	//	addBoolean("coexpression", false, "whether to calc coexperssion matrices (large memory for small bin size)");
	addBoolean("overwrite", true, "whether to overwrite existing results files");
		addBoolean("keepAlignment", false, "whether to keep alignment for MSA");
	//	addBoolean("attempt5rescue", false, "whether to attempt rescue of leader sequence if extra unmapped 5 read");
	//	addBoolean("attempt3rescue", false, "whether to attempt rescue of leader sequence if extra unmapped 5 read");
		//addBoolean("writePolyA", false, "whether write reads with polyA in middle");
		addBoolean("coronavirus", false, "whether to run in coronavirus mode (necessary to do breakpoint analysis, but takes more memory)");
		addBoolean("writeIsoforms", false, "whether to write isoforms");
		addBoolean("enforceKnownLinkedExons", false, "whether to use exon_exon boundaries from annotation ");
		addBoolean("writeH5", false,"whether to write h5 outputs");
		addString("mm2_path", "/sw/minimap2/current/minimap2",  "minimap2 path", false);
		addString("mm2Preset", "splice",  "preset for minimap2", false);
		addString("mm2_splicing", null, "splicing option", false);
		addInt("supplementaryQ", 30,"quality threshold for including supplementary mapping");
		addString("output",null, "Output file, default is stdout");
	//	addBoolean("writeBed", false, "whether to write bed",false);
		//addString("mm2Preset", "map-ont",  "preset for minimap2", false);
		addString("mm2_memory", (Runtime.getRuntime().maxMemory()-1000000000)+"",  "minimap2 memory", false);
		addInt("mm2_threads", 4, "threads for mm2", false);
		addStdHelp();
	}

	
 static void printParams(File resDir, String[] args1) throws IOException{
		File params_file = new File(resDir, "params.txt");
		PrintWriter pw = new PrintWriter(new FileWriter(params_file));
		Date d = new Date();
		
		 // System.getProperties().list(System.out);
		String mem = "Xmx="+(int)Math.ceil((double)Runtime.getRuntime().maxMemory()/(1e6))+"m";
		pw.println("User directory\n"+System.getProperty("user.dir"));
		pw.println("Date\n"+d.toString());
		pw.println("OS="+System.getProperty("os.name")+"\nJVM="+System.getProperty("java.version"));
		pw.println(mem);
		pw.println("Run commands");
		for(int i=0; i<args1.length; i+=1){
			pw.print(args1[i]+" ");
		}
		pw.close();
		
	}
 public static int supplementaryQ = 20;
 public static boolean coronavirus = true;
 public static int maxReads = Integer.MAX_VALUE;
 public static String pool_sep="";
 public static boolean limit_to_read_list = true;
	public static boolean sequential = true;
	public static boolean RNA;
	public static String mm2_index;
	public static String sampleID;
 public static void run(CommandLine cmdLine, String bamFiles,
		// String barcode_files,
		 String chrs, String chrsToIgnore,  boolean fastq, String reference) throws IOException{
		int qual = cmdLine.getIntVal("qual");
		//int bin = cmdLine.getIntVal("bin");
		
		CigarCluster.singleGFF = cmdLine.getBooleanVal("singleGFF");
		verbose=cmdLine.getBooleanVal("verbose");
		readsToSkip = cmdLine.getIntVal("readsToSkip");
		Outputs.overwrite  = cmdLine.getBooleanVal("overwrite");
		String coords_msa = cmdLine.getStringVal("coords_msa");
		if(coords_msa!=null && (new File(coords_msa)).exists()){
			BufferedReader br = new BufferedReader(new FileReader(coords_msa));
			String st = br.readLine();
			while((st = br.readLine())!=null){
				String[] coords = st.split(",");
				int[] pos = new int[2];
				//String[] coords1 = coords.split("_");
				pos[0] = Integer.parseInt(coords[1]);
				pos[1] = Integer.parseInt(coords[2]);
				IdentityProfile1.coords.add( pos);
				IdentityProfile1.coords_name.add(coords[0].replace(" ", "_"));
			
			}
			br.close();
		}
		Outputs.readsOutputFile = cmdLine.getStringVal("readsOutputFile");
		int breakThresh = cmdLine.getIntVal("breakThresh");
		String pattern = cmdLine.getStringVal("pattern");
		String readList = cmdLine.getStringVal("readList");
		//String genesToInclude = cmdLine.getStringVal("genesToInclude");
		
		//boolean overwrite  = cmdLine.getBooleanVal("overwrite");
		int startThresh = cmdLine.getIntVal("startThresh");
		int endThresh = cmdLine.getIntVal("endThresh");
		 maxReads = cmdLine.getIntVal("maxReads");
		//SequenceUtils.max_per_file = maxReads*4;
		 ViralTranscriptAnalysisCmd2.supplementaryQ =cmdLine.getIntVal("supplementaryQ");
		reverseForward = cmdLine.getBooleanVal("reverseForward");
		max_umi = cmdLine.getIntVal("max_umi");
		min_umi = cmdLine.getIntVal("min_umi");
		illumina=cmdLine.getBooleanVal("illumina");
		Outputs.url = cmdLine.getStringVal("api_url");
		String output = cmdLine.getStringVal("output", null);
		String output_join = cmdLine.getStringVal("join_out", null);
		String output_overlap = cmdLine.getStringVal("overlap_out", null);
		Outputs.format="json";
		if(output==null) {
			Outputs.outputstream = System.out;
		}else {
		  Outputs.outputstream =  output.endsWith(".gz")  ? new PrintStream(new GZIPOutputStream(new FileOutputStream(output))) : new PrintStream(new FileOutputStream(output));
		  Outputs.format = output.endsWith(".json.gz") || output.endsWith(".json") ? "json" : "tsv";
		}
		Outputs.joinOut =output_join==null ? null : (output_join.endsWith(".gz") ? new PrintStream(new GZIPOutputStream(new FileOutputStream(output_join))) : new PrintStream(new FileOutputStream(output_join)));
		Outputs.overlapOut =output_overlap==null ? null : (output_overlap.endsWith(".gz") ? new PrintStream(new GZIPOutputStream(new FileOutputStream(output_overlap))) : new PrintStream(new FileOutputStream(output_overlap)));

		Outputs.annotation_mode = Arrays.asList(cmdLine.getStringVal("annotation_mode").split(":"));
		Outputs.writeH5 = cmdLine.getBooleanVal("writeH5");
		Outputs.report = cmdLine.getIntVal("report");
		String[] round = cmdLine.getStringVal("round", "10").split(":");
		CigarHash2.round = new Integer[ round.length];
		for(int i=0; i<round.length;i++)CigarHash2.round[i] = Integer.parseInt(round[i]);
	//	CigarHash2.round1 = cmdLine.getIntVal("bin1");
		String types_to_include = cmdLine.getStringVal("types_to_include", null);
		if(types_to_include!=null) IdentityProfile1.types_to_include = Arrays.asList(types_to_include.split(":"));
		Annotation.enforceStrand = cmdLine.getBooleanVal("enforceStrand");
		IdentityProfile1.trainStrand = cmdLine.getBooleanVal("trainStrand");
		GFFAnnotation.enforceKnownLinkedExons = cmdLine.getBooleanVal("enforceKnownLinkedExons");
		if(IdentityProfile1.trainStrand && Annotation.enforceStrand) throw new RuntimeException("cant enforce and train");
		String RNAstr = 		cmdLine.getStringVal("RNA","true");
		//RNA = new boolean[bamFiles.length];
		Outputs.sampleName =bamFiles;
		if(RNAstr!=null){
				if(RNAstr.equals("name")) {
						RNA = bamFiles.toLowerCase().indexOf("rna")>0;
				}
				else RNA= Boolean.parseBoolean(RNAstr) ;
		}else{
			RNA=true;
		}
		System.err.println(Arrays.asList(RNA));
	//	if(IdentityProfile1.trainStrand &&
	//	Outputs.executor=  
	//			cmdLine.getIntVal("max_threadsIO")==1 ? 
	//			Executors.newSingleThreadExecutor():  Executors.newFixedThreadPool(cmdLine.getIntVal("maxThreads"));
		IdentityProfileHolder.executor=  cmdLine.getIntVal("maxThreads")==1 ? null:  Executors.newFixedThreadPool(cmdLine.getIntVal("maxThreads"));
		IdentityProfile1.qual_thresh = cmdLine.getDoubleVal("leftOverQualThresh");
		//ViralTranscriptAnalysisCmd2.combineOutput = cmdLine.getBooleanVal("combineOutput");
		String d_thresh = cmdLine.getStringVal("isoformDepthThresh");
		int isoformDepthThresh  = Integer.parseInt(d_thresh);
		int coverageDepthThresh = cmdLine.getIntVal("coverageDepthThresh");
		IdentityProfile1.msaDepthThresh =(int) Math.floor(cmdLine.getDoubleVal("msaDepthThresh"));
	//IdentityProfile1.extra_threshold = cmdLine.getIntVal("extra_threshold");
//	IdentityProfile1.tryComplementOnExtra = cmdLine.getBooleanVal("tryComplementOnExtra");
//	IdentityProfile1.reAlignExtra = cmdLine.getBooleanVal("reAlignExtra");
//	IdentityProfile1.attempt5rescue = cmdLine.getBooleanVal("attempt5rescue");
//	IdentityProfile1.attempt3rescue = cmdLine.getBooleanVal("attempt3rescue");
	String[]	fail_thresh = cmdLine.getStringVal("fail_thresh").split(":");
	ViralTranscriptAnalysisCmd2.fail_thresh = Double.parseDouble(fail_thresh[0]);
	ViralTranscriptAnalysisCmd2.fail_thresh1 = fail_thresh.length>1 ? Double.parseDouble(fail_thresh[1]) : Double.parseDouble(fail_thresh[0]);

		Pattern patt = Pattern.compile(":");
		IdentityProfile1.includeStartEnd = cmdLine.getBooleanVal("includeStart");
	SequenceOutputStream1.max_seqs_per_cluster = cmdLine.getIntVal("max_seqs_per_cluster");
		 coronavirus = cmdLine.getBooleanVal("coronavirus");
		String[] msaOpts = cmdLine.getStringVal("doMSA").split(":"); //e.g 5_3:sep or all:sep
		//String msa_source = cmdLine.getStringVal("msa_source");
		//double probInclude = cmdLine.getDoubleVal("probInclude");
		/*if(msa_source!=null){
			if(msa_source.equals("all")){
				for(int j=0; j<bamFiles.length;j++){
					Outputs.msa_sources.put(j, j);
				}
			}else{
				String[] str = msa_source.split(";");
				outer: for(int j=0; j<bamFiles.length;j++){
					for(int i=0; i<str.length;i++){
							String[] str1 = str[i].split(",");
							for(int k=0; k<str1.length; k++){
								if(bamFiles[j].indexOf(str1[k])>=0){
									Outputs.msa_sources.put(j, i);
									continue outer;
								}
							}
						}
				}
			}
			System.err.println("msa_sources "+Outputs.msa_sources);
		}else{
			
			for(int i=0; i<len; i++){
				Outputs.msa_sources.put(i, 0);
			}
		}
		for(int i=0; i<msaOpts.length; i++){
			if(msaOpts[i].startsWith("span=")) msaOpts[i] = msaOpts[i].split("=")[1];
		}*/
		Outputs.doMSA =Arrays.asList((msaOpts[0].equals("all") ?"5_3,no5_3,5_no3,no5_no3": msaOpts[0]).split(",")) ;
		Outputs.minClusterEntries = cmdLine.getIntVal("minClusterEntries");
		if(msaOpts[0].equals("false")){
			Outputs.doMSA = null;
		}
		Outputs.keepAlignment = cmdLine.getBooleanVal("keepAlignment");
		Outputs.keepinputFasta = Outputs.keepAlignment ;
		String chromsToRemap = cmdLine.getStringVal("chromsToRemap");
		
		//Outputs.MSA_at_cluster = false;
		Outputs.calcBreaks=cmdLine.getBooleanVal("calcBreaks");// whether to calculate data for the break point heatmap, true for SARS_COV2
	//	boolean filterBy5_3 = false;// should be true for SARS_COV2
		boolean annotByBreakPosition = false;  // should be true for SARS_COV2
		//Outputs.writePolyA = cmdLine.getBooleanVal("writePolyA");
		CigarCluster.recordDepthByPosition = cmdLine.getBooleanVal("recordDepthByPosition");
		CigarCluster.recordStartEnd = cmdLine.getBooleanVal("recordDepthByPosition");
		IdentityProfileHolder.overlap_max = cmdLine.getIntVal("overlap_max");
		IdentityProfileHolder.gap_max = cmdLine.getIntVal("gap_max");
		IdentityProfileHolder.overlap_thresh = cmdLine.getDoubleVal("overlap_thresh");
		IdentityProfile1.min_first_last_exon_length = cmdLine.getIntVal("min_first_last_exon");
		if(coronavirus){
			System.err.println("running in coronavirus mode");
			SequenceUtils.mm2_splicing= "-"+cmdLine.getStringVal("mm2_splicing","un");
		//sequential = true;
		}else{
		//	TranscriptUtils.coronavirus = false;
			SequenceUtils.mm2_splicing= "-"+cmdLine.getStringVal("mm2_splicing","uf");
		}
	//	SequenceUtils.mm2Preset=fnull;//"map-ont";
		sequential = cmdLine.getBooleanVal("sequential");
		sorted = cmdLine.getBooleanVal("sorted");
		

		Outputs.library = new File(cmdLine.getStringVal("library"));
		
		String annotFile = cmdLine.getStringVal("annotation_file");
		/*	static void errorAnalysis(String bamFile_,
		 String refFile,  String readList,    String resdir, 
		String pattern, int qual, int round, 
		int break_thresh, int startThresh, int endThresh, int max_reads, 
		boolean calcBreaks1 ,  boolean annotByBreakPosition, String chrToInclude, String chrToIgnore, 
		int writeIsoformDepthThresh, int writeCoverageDepthThresh, boolean fastq, String chromsToRemap, 
		String annot_file) throws IOException {*/

			errorAnalysis(bamFiles,reference, readList,
				pattern, qual,  breakThresh, startThresh, endThresh,maxReads,  
				Outputs.calcBreaks, annotByBreakPosition,  chrs, chrsToIgnore,  isoformDepthThresh, coverageDepthThresh,  fastq, 
				chromsToRemap, annotFile);
	}
 public static boolean sorted = true;
public static boolean allowSuppAlignments = true;; // this has to be true for allowing supp alignments

 static double tme0;
 public static String[] readOpts(CommandLine cmdLine, String[] args1, String  opts_file) throws IOException{
	 String[] args = cmdLine.stdParseLine(args1);
	 Map<String, String> keyval = new HashMap<String, String>();
		String optsFiles  = cmdLine.getStringVal(opts_file);
		String optsType = cmdLine.getStringVal("optsType");
		System.err.println("optsType "+optsType);
		
		if(optsFiles!=null){
			BufferedReader br = new BufferedReader(new FileReader(optsFiles));
			String st = "";
			int column = 0;
			while((st=br.readLine())!=null){
				if(st.trim().startsWith("#") || st.length()==0) continue;
				String st1 = st.trim().split("#")[0].trim();
				String[] st2 = st1.split(";");
				if(st1.contains("optsType")){
					if(optsType!=null){
						column = -1;
						inner: for(int j=0; j<st2.length; j++){
							String[] st3 = st2[j].trim().split("=");
							if(st3[1].contains(optsType)){
								System.err.println("using column "+j+" where specified");
								column=j;
								break inner;
							}
						}
						if(column <0){
							throw new RuntimeException("!! no column "+optsType);
						}
					}
					continue;
				}
				{
					int k  = column < st2.length ? column : 0; 
					String[] st3 = st2[k].trim().split("=");
					if(st3.length!=2) {
						System.err.println(st1);
						System.err.println(st2[k].trim());
						throw new RuntimeException("wrong format, needs to be key=val "+st2[k]);
					}
					if(!keyval.containsKey(st3[0])) {
						keyval.put(st3[0], st3[1]);
					}
				}
			}
			br.close();
			
		}
	 // command line parameters last to override
	 for(int i=0; i<args1.length; i++){
			String[] st2 = args1[i].split("=");
			if(st2.length<2) throw new RuntimeException("wrong format, needs to be ke=val "+args1[i]);
			keyval.put(st2[0], st2[1]);
		}
	 args = new String[keyval.size()];
		Iterator<Entry<String, String>> it = keyval.entrySet().iterator();
		for(int i=0; it.hasNext(); i++){
			Entry<String, String> ent = it.next();
			args[i] = ent.getKey()+"="+ent.getValue();
		}
		
	 return cmdLine.stdParseLine(args);
 }
	public static void main(String[] args1) throws IOException, InterruptedException {
		tme0 = System.currentTimeMillis();
		
		CommandLine cmdLine = new ViralTranscriptAnalysisCmd2();
		String[] args2 = ViralTranscriptAnalysisCmd2.readOpts(cmdLine, args1,"optsFile");
		
		//String fastqFile = cmdLine.getStringVal("fastqFile");
		
		
		//Outputs.writeBed  = cmdLine.getBooleanVal("writeBed");
		Outputs.writeIsoforms = cmdLine.getBooleanVal("writeIsoforms");
barcode_list = cmdLine.getStringVal("barcode_list");
barcode_file = cmdLine.getStringVal("barcode_file");
		SequenceUtils.mm2_threads = cmdLine.getIntVal("mm2_threads");
		SequenceUtils.mm2_mem = cmdLine.getStringVal("mm2_mem");
		SequenceUtils.mm2_path = cmdLine.getStringVal("mm2_path");
		SequenceUtils.mm2Preset = cmdLine.getStringVal("mm2Preset");
		
		 Barcodes.barcode_extent = cmdLine.getIntVal("barcode_extent");
		 Barcodes.barcode_ignore = cmdLine.getIntVal("barcode_ignore");
		 Barcodes.tolerance_barcode = cmdLine.getIntVal("tolerance_barcode");
		// ViralTranscriptAnalysisCmd2.barcode_list = cmdLine.getStringVal("barcode_list");
		 
		// PolyAT.set_bc_len(cmdLine.getIntVal("bc_len_AT"));
		// PolyAT.edit_thresh_AT = cmdLine.getIntVal("edit_thresh_AT");
		// PolyAT.max_dist_to_end_AT = cmdLine.getIntVal("max_dist_to_end_AT");
		// PolyAT.dist_to_ignore_AT = cmdLine.getIntVal("dist_to_ignore_AT");


		ViralTranscriptAnalysisCmd2.exclude_indeterminate_strand= cmdLine.getBooleanVal("exclude_indeterminate_strand");
		ViralTranscriptAnalysisCmd2.exclude_polyT_strand= cmdLine.getBooleanVal("exclude_polyT_strand");
		ViralTranscriptAnalysisCmd2.exclude_reads_without_barcode= cmdLine.getBooleanVal("exclude_reads_without_barcode");

		
		String reference = cmdLine.getStringVal("reference");
		File refFile = reference==null ? null : new File(reference);
	//	if(!refFile.exists()) throw new RuntimeException("ref file does not exist");
		
		
		
		
	//	String resdir = cmdLine.getStringVal("resdir");
		
		
		String chroms= cmdLine.getStringVal("chroms_to_include");
		String chroms_ignore= cmdLine.getStringVal("chroms_to_ignore");
		
		//boolean gff = annot_file !=null && (annot_file.contains(".gff") || annot_file.contains(".gtf"));
		
		//File resDir = new File(resdir);
		//if(!resDir.exists())resDir.mkdir();
		//printParams(resDir, args1);
		String inputFiles_ = cmdLine.getStringVal("inputFile","-");
		//if(inputFile.equals("-"))stdin=true;
		String[] pathToInpFile = inputFiles_.split("/");
		sampleID = cmdLine.getStringVal("sampleID",pathToInpFile[pathToInpFile.length-1]);
			boolean bam = true; //inputFiles_.equals("-") ||  inputFiles_.endsWith(".bam") || inputFiles_.endsWith("sam") ;
			//String[] barcode_files = new String[inputFiles_.length];
			/*if(barcode_file!=null){
				barcode_files = barcode_file.split(":");
				File fi = new File(barcode_files[0]);
				if(!fi.exists()) {
					throw new RuntimeException(" barcode file does not exist "+fi.getAbsolutePath());
				}
				if(barcode_files.length!=inputFiles_.length){
					throw new RuntimeException(" wrong lenght of barcodes");
				}
			}else if(barcode_list!=null){
			File barcode_file_ = new File(barcode_list);
			if(barcode_file_.exists()){
				if(inputFiles_[0].equals("-")){
					barcode_files[0] = barcode_list;
				}else{
				BufferedReader br = new BufferedReader(new FileReader(barcode_file_));
				String st = "";
				List<String> ifs = new ArrayList<String>();
				while((st=br.readLine())!=null){
					if(st.length()>0 && !st.trim().startsWith("#")){
						String[] str = st.trim().split("\t");
						for(int i=0; i<inputFiles_.length; i++){
							File fi = new File(inputFiles_[i]);
							if(fi.getName().equals(str[0])){
								barcode_files[i] = str[1];
							}
						}
						
					}
				}
				}
			}
			}*/
			
			
			boolean fastq = !bam;
			if(fastq){
				try{
					//make a minimap index
					boolean saveSeqs = true;
				mm2_index = SequenceUtils.minimapIndex(refFile,  false, saveSeqs);
				}catch(Exception exc){
					exc.printStackTrace();
				}
			}
		run(cmdLine, inputFiles_,     chroms, chroms_ignore, fastq, reference);
		System.err.println("shutting down executors");
		
		IdentityProfileHolder.shutDownExecutor();
		Outputs.shutDownExecutor();
		double time1 = System.currentTimeMillis();
		System.err.println((time1-tme0)/(1000)+ " seconds");
		// paramEst(bamFile, reference, qual);
	}
//public static boolean combineOutput = false;
	public static double fail_thresh = 7.0;
	public static double fail_thresh1 = 14.0;
	public static String barcode_list=null;
	public static String barcode_file=null;
	public static Barcodes[] barcodes  = null;
	/**
	 * Error analysis of a bam file. Assume it has been sorted
	 */
	static void errorAnalysis(String bamFile_,
			 String refFile,  String readList,  
			String pattern, int qual, 
			int break_thresh, int startThresh, int endThresh, int max_reads, 
			boolean calcBreaks1 ,  boolean annotByBreakPosition, String chrToInclude, String chrToIgnore, 
			int writeIsoformDepthThresh, int writeCoverageDepthThresh, boolean fastq, String chromsToRemap, 
			String annot_file) throws IOException {
		boolean cluster_reads = true;
	
	//	Annotation.tolerance = round;
		
		IdentityProfile1.annotByBreakPosition = annotByBreakPosition;
		CigarHash.cluster_by_annotation =true;// cluster_by_annotation;
		TranscriptUtils.startThresh = startThresh;
		TranscriptUtils.endThresh = endThresh;
		IdentityProfile1.writeIsoformDepthThresh =writeIsoformDepthThresh;
		IdentityProfile1.writeCoverageDepthThresh =writeCoverageDepthThresh;
		Map<String, Integer> reads= null;//new HashMap<String, Integer>();
	
		IdentityProfile1.break_thresh = break_thresh;
		System.err.println("reading genomes " + (System.currentTimeMillis()- tme0)/1000+" secs");
		ArrayList<Sequence> genomes = refFile==null ? null : SequenceReader.readAll(refFile, Alphabet.DNA());
		System.err.println("done reading genomes" + (System.currentTimeMillis()- tme0)/1000+" secs");
		Map<String, Integer>chrNameToIndex = genomes==null ? null : new HashMap<String, Integer>();
		if(genomes!=null){
			for(int i=0; i<genomes.size(); i++){
					chrNameToIndex.put(genomes.get(i).getName(),i);
			}
		}
		
	
		String in_nmes;
			if(!bamFile_.equals("-")){
				File bam = new File( bamFile_);
				in_nmes = bam.getName().split("\\.")[0];
			}else{
				in_nmes = "stdin";
			}
		SamReader samReaders= null;// = new SamReader[len];
		boolean allNull = true;
		Iterator<SAMRecord> samIter=null;
			SamReaderFactory.setDefaultValidationStringency(ValidationStringency.SILENT);
			if(bamFile_.equals("-")){
				System.err.println("reading from stdin");
				samReaders = SamReaderFactory.makeDefault().open(SamInputResource.of(System.in));
			}else{

				samReaders = SamReaderFactory.makeDefault().open(new File(bamFile_));
			}
			Map<String, Object> m = new HashMap();
			Map<String, Object> flags = new HashMap();
			List<SAMProgramRecord> pr = samReaders.getFileHeader().getProgramRecords();
			flags.put("program_records",pr );
			List<String> programs = new ArrayList<String>();
			for(int kk=0; kk<pr.size(); kk++) {
				SAMProgramRecord pr1 = pr.get(kk);
				programs.add(pr1.getProgramName());
				if(pr1.getCommandLine().startsWith("minimap2")) {
					List<String>str = Arrays.asList(pr1.getCommandLine().split("\\s+"));
					StringBuffer sb = new StringBuffer();
					for(int jk=1; jk<str.size()-2; jk++) sb.append(str.get(jk)+" ");
					flags.put("minimap2_options", sb.toString().trim());
					flags.put("minimap2_version", pr1.getAttribute("VN"));
					flags.put("minimap2_reference", str.get(str.size()-2).replaceAll(".mmi", ""));
					flags.put("fastq", str.get(str.size()-1));
					if(sampleID==null || sampleID.equals("-") ||sampleID.equals("fastq")) {
						sampleID = (String) flags.get("fastq");
					}
				}
			}
			flags.put("programs", programs);
			flags.put("annot", Outputs.annotation_mode);
			flags.put("round", CigarHash2.round);
			flags.put("br", IdentityProfile1.break_thresh);
			m.put("id", sampleID);
			m.put("flags", flags);
				samIter = samReaders.iterator();
				Gson gson = new Gson();
				String json_str = gson.toJson(m);
				Outputs.outputstream.print(json_str);Outputs.outputstream.println();
			allNull = false;

			Map<String, Object> expt=new HashMap<String, Object>() ;
			Map<String, Object> flags1 = new HashMap<String, Object>();
			expt.put("sampleID",sampleID); // prob not best sampleID
			expt.put("flags",flags1);
			flags1.put("RNA", ViralTranscriptAnalysisCmd2.RNA);
			Outputs 	outp = new Outputs( in_nmes,   true, CigarCluster.recordDepthByPosition,expt); 
			
			
			
			Annotation 	annot = annot_file == null ? new EmptyAnnotation( ) : 
					new Annotation(new File(annot_file), 1);
			
			
			IdentityProfileHolder profile =  	
					new IdentityProfileHolder(genomes, null, outp,  in_nmes, calcBreaks1, annot);
			
			
			int numReads = 0;

			int numNotAligned = 0;
			int numSecondary =0;
			
			long time0 = System.currentTimeMillis();//- tme0)/1000.0
			int cnt0=0;
			int ij=0;
			String previousRead = "";
			int previousRefIndex=-1;
			Set<String> skipped = new HashSet<String>();
			List<SAMRecord> sams = new ArrayList<SAMRecord>();
			int primaryIndex=-1;
			String prev_readnme = "";
			System.err.println("RNA "+Arrays.asList(RNA));
		
			outer: for ( ij=0; samIter.hasNext() ;ij++ ) {
				
				
				final SAMRecord sam=samIter.next();
				
			    if(sam==null){
			    	break outer;
			    }
			    if(outp.skip(sam.getReadName())){
					continue outer;
				}
			    if(verbose){
					System.err.println(sam.getReadName());
				}
				if (sam.getReadUnmappedFlag()) {
				//	outp.writeUnmapped(sam.getReadName());
					numNotAligned++;
					continue;
				}
				
				byte[] b = sam.getBaseQualities();
				double q1 = 500;
				if(b.length>0){
					double sump = 0;
					for(int i=0; i<b.length; i++){
						sump+= Math.pow(10, -b[i]/10.0);
					}
					sump = sump/(double) b.length;
					 q1 = -10.0*Math.log10(sump);
					
					if(q1 < fail_thresh) {
						//System.err.println(q1);
						continue;
					}
				}
				
				if(sam==null) break outer;

				numReads++;
				if (sam.getMappingQuality() < qual) {
					numNotAligned++;
					continue;
				}
				boolean supplementary = false;
				if(sam.isSecondaryOrSupplementary()) {
					if(sam.isSecondaryAlignment()) continue;
					if(!ViralTranscriptAnalysisCmd2.allowSuppAlignments) continue;
					if( sam.getReadName().equals(previousRead) ){
						supplementary = true;
						int qual1 = sam.getMappingQuality();
						
						if(qual1 < supplementaryQ  ){
							continue; // need secondary read quality of 5
						}else{
							if(sorted) throw new RuntimeException("!!");;
						}
						
					}else{
					numSecondary++;
					
					continue;
					}
				}
				
				
					if(sams.size()>0 && (!sam.getReadName().equals(previousRead) )){
						try{
							profile.identity1(new ArrayList<SAMRecord>(sams),  cluster_reads,     genomes, primaryIndex);
							sams.clear();
							primaryIndex=-1;
						}catch(NumberFormatException exc){
						
							exc.printStackTrace();
						}
					}
					if(!sam.isSecondaryOrSupplementary()){
						primaryIndex = sams.size();
						
					}
					sams.add(sam);

					previousRead = sam.getReadName();
					previousRefIndex = sam.getReferenceIndex();
			}
			if(sams.size()>0){
				try{
					
					profile.identity1(new ArrayList<SAMRecord>(sams),  cluster_reads,     genomes, primaryIndex);
					sams.clear();
				}catch(NumberFormatException exc){
					exc.printStackTrace();
				}
			}
			if(samReaders!=null)samReaders.close();
			
		
		
			if(profile!=null){
				IdentityProfileHolder.waitOnThreads(300);
			}
			if(outp!=null) outp.close();
	
			
				
			{
				long delta =  System.currentTimeMillis() -  time0;
				double timeperread =  ((double)delta)/ ((double)ij);
				System.err.println(timeperread+" milliseconds per read at "+ij);
			}	
	}

	
}
 
