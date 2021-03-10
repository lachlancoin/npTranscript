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
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.FilenameFilter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Date;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.Executors;
import java.util.regex.Pattern;
import java.util.stream.Collectors;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;
import java.util.zip.ZipFile;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.fastq.BasicFastqWriter;
import htsjdk.samtools.fastq.FastqRecord;
import htsjdk.samtools.fastq.FastqWriter;
import htsjdk.samtools.util.SequenceUtil;
import japsa.seq.Alphabet;
import japsa.seq.Sequence;
import japsa.seq.SequenceReader;
import japsa.seq.ZipGFF;
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
import npTranscript.cluster.CigarCluster.Count;

/**
 * @author Lachlan Coin
 *
 */

@Deployable(scriptName = "npTranscript.run", scriptDesc = "Analysis of coronavirus sequence data")
public class ViralTranscriptAnalysisCmd2 extends CommandLine {
	

//	private static final Logger LOG = LoggerFactory.getLogger(HTSErrorAnalysisCmd.class);

public static String getAnnotationsToInclude(String annotationType, boolean useExons){
	if(annotationType!=null) return annotationType;//.split(":");
	 if(!useExons) {
		 return  "gene:ncRNA:pseudogene:miRNA";	
	 }
	 else{
		 return "all";
//		 return "exon";
	 }
		
}
	public ViralTranscriptAnalysisCmd2() {
		super();
		Deployable annotation = getClass().getAnnotation(Deployable.class);
		setUsage(annotation.scriptName() + " [options]");
		setDesc(annotation.scriptDesc());
		addString("bamFile", null, "Name of bam file", false);
		addString("fastqFile", null, "Name of bam file", false);
		addString("reference", null, "Name of reference genome", true);
		addString("annotation", null, "ORF annotation file or GFF file", false);
		addBoolean("useExons", true, "wehether to use exons");
		addBoolean("baseBreakPointsOnFirst", false, "whether the break points should be preferentially based on the first bam");
		addBoolean("sorted", false, "whether bamfile sorted in terms of position.");
		addBoolean("sequential", true, "whether to iterate through one bam at a time (in sequence) or in parallel.");
		addBoolean("firstIsTranscriptome", false, "whether first bam is transcriptome");
		addBoolean("enforceStrand", false, "whether to enforce strandedness");
		addString("readList", "", "List of reads", false);
		addString("gffThresh","1", "reports if greater than in dataset");

		addInt("maxTranscriptsPerGeneInGFF",1,"Maximum number of transcripts per gene (highest abundance first");
			addString("annotType", null, "Type of annotation (only included if annotation is GFF file", false);
		addString("chroms_to_include", "all", "Restrict to these chroms, colon delimited", false);
		addString("chroms_to_ignore", "none", "Ignore these chroms", false);
	//	addString("bedChr", null, "Use this for the chrom in bed chr, e.g. NC_045512v2, false");
        addBoolean("singleGFF", false, "whether to have single output gff, rather than different for different input. ");
		addString("resdir", "results"+System.currentTimeMillis(), "results directory");
		addString("GFF_features", "gene_name:description:gene_id:gene_biotype:gene_id", "GFF feature names");
		addString("RNA", "name", "If is direct RNA.  Can be tab delimmited boolean, e.g. true:true  or if set to name it will look for RNA string in name");
		addBoolean("trainStrand", false,"whether to produce training data for strand correction");
		addBoolean("recordDepthByPosition", false, "whether to store position specific depth (high memory");
	//	addString("annotToRemoveFromGFF",null, "annotation to remove from GFF , BED and ref files");
		addInt("maxReads", Integer.MAX_VALUE, "ORF annotation file");
		addDouble("probInclude", 1.0, "probability of including each read");
		addInt("minClusterEntries",10,"threshold for consensus");
	//	addBoolean("tryComplementOnExtra", false, "look for negative strand matches on left over seqs");
		addBoolean("reAlignExtra", false, "whether to try realigning the extra sequence");
	//	addBoolean("combineOutput", false, "whether to combine output from different chroms");
		addString("pattern", null, "Pattern of read name, used for filtering");
		addString("span", "protein_coding", "Filtering span.  Use all not to filter.");
		addInt("qual", 0, "Minimum quality required");
		addInt("bin", 1, "Bin size for numerical hashing");
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
		//addString("aligner", "kalign" , "Options: kalign3, poa, spoa, abpoa");
		addInt("startThresh", 100, "Threshold for having 5'");
		addInt("endThresh", 100, "Threshold for having 3'");
		addInt("maxThreads", 8, "max threads for processing");
		addInt("extra_threshold", 200, "Threshold saving umatched 3'or 5'parts of reads");
		//addDouble("overlapThresh", 0.95, "Threshold for overlapping clusters");
	//	addBoolean("coexpression", false, "whether to calc coexperssion matrices (large memory for small bin size)");
	//	addBoolean("overwrite", false, "whether to overwrite existing results files");
		addBoolean("keepAlignment", false, "whether to keep alignment for MSA");
		addBoolean("attempt5rescue", false, "whether to attempt rescue of leader sequence if extra unmapped 5 read");
		addBoolean("attempt3rescue", false, "whether to attempt rescue of leader sequence if extra unmapped 5 read");
		addBoolean("writePolyA", false, "whether write reads with polyA in middle");
		addBoolean("coronavirus", true, "whether to run in coronavirus mode (necessary to do breakpoint analysis, but takes more memory)");
		addBoolean("writeGFF", false, "whether to output gff ");
		addBoolean("writeIsoforms", false, "whether to write isoforms");
		addBoolean("enforceKnownLinkedExons", false, "whether to use exon_exon boundaries from annotation ");
		addBoolean("writeH5", false,"whether to write h5 outputs");
		addString("mm2_path", "/sw/minimap2/current/minimap2",  "minimap2 path", false);
		addString("mm2Preset", "splice",  "preset for minimap2", false);
		addInt("supplementaryQ", 1000,"quality threshold for including supplementary mapping");
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
 public static boolean coronavirus = false;
 public static int maxReads = Integer.MAX_VALUE;
 public static String pool_sep="";
 public static boolean limit_to_read_list = true;
	public static boolean sequential = false;
	public static boolean[] RNA;
	public static String mm2_index;
 public static void run(CommandLine cmdLine, String[] bamFiles, String resDir,File anno, String chrs, String chrsToIgnore,  boolean fastq, String reference) throws IOException{
		Count.baseBreakPointsOnFirst = cmdLine.getBooleanVal("baseBreakPointsOnFirst");
		int qual = cmdLine.getIntVal("qual");
		int bin = cmdLine.getIntVal("bin");
		CigarCluster.singleGFF = cmdLine.getBooleanVal("singleGFF");
		int breakThresh = cmdLine.getIntVal("breakThresh");
		String pattern = cmdLine.getStringVal("pattern");
		String annotFile = cmdLine.getStringVal("annotation");
		String[] readList = cmdLine.getStringVal("readList").split(":");
		//String genesToInclude = cmdLine.getStringVal("genesToInclude");
		String  annotationType = getAnnotationsToInclude(cmdLine.getStringVal("annotType"), cmdLine.getBooleanVal("useExons"));
		//boolean overwrite  = cmdLine.getBooleanVal("overwrite");
		int startThresh = cmdLine.getIntVal("startThresh");
		int endThresh = cmdLine.getIntVal("endThresh");
		 maxReads = cmdLine.getIntVal("maxReads");
		//SequenceUtils.max_per_file = maxReads*4;
		 ViralTranscriptAnalysisCmd2.supplementaryQ =cmdLine.getIntVal("supplementaryQ");
		
		Outputs.maxTranscriptsPerGeneInGFF = cmdLine.getIntVal("maxTranscriptsPerGeneInGFF");
		Outputs.writeH5 = cmdLine.getBooleanVal("writeH5");
		Annotation.enforceStrand = cmdLine.getBooleanVal("enforceStrand");
		IdentityProfile1.trainStrand = cmdLine.getBooleanVal("trainStrand");
		GFFAnnotation.enforceKnownLinkedExons = cmdLine.getBooleanVal("enforceKnownLinkedExons");
		if(IdentityProfile1.trainStrand && Annotation.enforceStrand) throw new RuntimeException("cant enforce and train");
		String[] RNAstr = 		cmdLine.getStringVal("RNA").split(":");
		RNA = new boolean[bamFiles.length];
		if(RNAstr!=null){
			if(RNAstr.length==1){
				if(RNAstr[0].equals("name")) {
					for(int i=0; i<bamFiles.length; i++){
						RNA[i] = bamFiles[i].indexOf("RNA")>0;
					}
				}
				else Arrays.fill(RNA, Boolean.parseBoolean(RNAstr[0]) );
			}
			else for(int i=0; i<RNAstr.length; i++) RNA[i] = Boolean.parseBoolean(RNAstr[i]);
		}else{
			Arrays.fill(RNA, false);
		}
	//	if(IdentityProfile1.trainStrand &&
	//	Outputs.executor=  
	//			cmdLine.getIntVal("max_threadsIO")==1 ? 
	//			Executors.newSingleThreadExecutor():  Executors.newFixedThreadPool(cmdLine.getIntVal("maxThreads"));
		IdentityProfileHolder.executor=  cmdLine.getIntVal("maxThreads")==1 ? null:  Executors.newFixedThreadPool(cmdLine.getIntVal("maxThreads"));
		IdentityProfile1.qual_thresh = cmdLine.getDoubleVal("leftOverQualThresh");
		//ViralTranscriptAnalysisCmd2.combineOutput = cmdLine.getBooleanVal("combineOutput");
		String[] d_thresh = cmdLine.getStringVal("isoformDepthThresh").split(":");
		int[] isoformDepthThresh  = new int[d_thresh.length];
		for(int i=0; i<d_thresh.length; i++){
			isoformDepthThresh[i] = Integer.parseInt(d_thresh[i]);
		}
		int coverageDepthThresh = cmdLine.getIntVal("coverageDepthThresh");
		IdentityProfile1.msaDepthThresh =(int) Math.floor(cmdLine.getDoubleVal("msaDepthThresh"));
	IdentityProfile1.extra_threshold = cmdLine.getIntVal("extra_threshold");
//	IdentityProfile1.tryComplementOnExtra = cmdLine.getBooleanVal("tryComplementOnExtra");
	IdentityProfile1.reAlignExtra = cmdLine.getBooleanVal("reAlignExtra");
	IdentityProfile1.attempt5rescue = cmdLine.getBooleanVal("attempt5rescue");
	IdentityProfile1.attempt3rescue = cmdLine.getBooleanVal("attempt3rescue");
	String[]	fail_thresh = cmdLine.getStringVal("fail_thresh").split(":");
	ViralTranscriptAnalysisCmd2.fail_thresh = Double.parseDouble(fail_thresh[0]);
	ViralTranscriptAnalysisCmd2.fail_thresh1 = fail_thresh.length>1 ? Double.parseDouble(fail_thresh[1]) : Double.parseDouble(fail_thresh[0]);

		Pattern patt = Pattern.compile(":");
		Outputs.numExonsMSA = cmdLine.getStringVal("numExonsMSA")=="none"  ?  Arrays.asList(new Integer[0]) : 
				patt.splitAsStream(cmdLine.getStringVal("numExonsMSA"))
		                            .map(Integer::valueOf)
		                            .collect(Collectors.toList());
		IdentityProfile1.includeStartEnd = cmdLine.getBooleanVal("includeStart");
		GFFAnnotation.setGFFFeatureNames(cmdLine.getStringVal("GFF_features").split(":"));
		GFFAnnotation.span_only = cmdLine.getStringVal("span").equals("all") ?new ArrayList<String>() :   Arrays.asList(cmdLine.getStringVal("span").split(":"));
		SequenceOutputStream1.max_seqs_per_cluster = cmdLine.getIntVal("max_seqs_per_cluster");
		 coronavirus = cmdLine.getBooleanVal("coronavirus");
		String[] msaOpts = cmdLine.getStringVal("doMSA").split(":"); //e.g 5_3:sep or all:sep
		String msa_source = cmdLine.getStringVal("msa_source");
		double probInclude = cmdLine.getDoubleVal("probInclude");
	//	String[] bamFiles =bamFile.split(":"); 
		int len =  bamFiles.length;
		if(msa_source!=null){
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
			System.err.println("msa_sources "+Outputs.msa_sources);
		}else{
			
			for(int i=0; i<len; i++){
				Outputs.msa_sources.put(i, 0);
			}
		}
		for(int i=0; i<msaOpts.length; i++){
			if(msaOpts[i].startsWith("span=")) msaOpts[i] = msaOpts[i].split("=")[1];
		}
		Outputs.doMSA =Arrays.asList((msaOpts[0].equals("all") ?"5_3,no5_3,5_no3,no5_no3": msaOpts[0]).split(",")) ;
		Outputs.minClusterEntries = cmdLine.getIntVal("minClusterEntries");
		if(msaOpts[0].equals("false")){
			Outputs.doMSA = null;
		}
		Outputs.keepAlignment = cmdLine.getBooleanVal("keepAlignment");
		Outputs.keepinputFasta = Outputs.keepAlignment ;
		String chromsToRemap = cmdLine.getStringVal("chromsToRemap");
		
		//Outputs.MSA_at_cluster = false;
		boolean calcBreaks=true;// = cmdLine.getBooleanVal("calcBreaks");// whether to calculate data for the break point heatmap, true for SARS_COV2
		boolean filterBy5_3 = false;// should be true for SARS_COV2
		boolean annotByBreakPosition = false;  // should be true for SARS_COV2
		Outputs.writePolyA = cmdLine.getBooleanVal("writePolyA");
		CigarCluster.recordDepthByPosition = cmdLine.getBooleanVal("recordDepthByPosition");
		CigarCluster.recordStartEnd = cmdLine.getBooleanVal("recordDepthByPosition");

		if(coronavirus){
		//	
		//	CigarCluster.recordDepthByPosition = true;
			IdentityProfile1.min_first_last_exon_length = 0;
			System.err.println("running in coronavirus mode");
			calcBreaks  = true; 
			filterBy5_3 = true;
		//	Outputs.writeGFF=true;
			Outputs.calcBreaks=true;
		//	Outputs.MSA_at_cluster = true;
			TranscriptUtils.checkAlign = true;
			TranscriptUtils.coronavirus = true;
			IdentityProfile1.extra_threshold1 =50;
			annotByBreakPosition = true;
			TranscriptUtils.writeAnnotP = true;
			
		//	CigarHash2.subclusterBasedOnStEnd = false;
			SequenceUtils.mm2_splicing= "-un";
		//sequential = true;
		}else{
			//sequential = false;
		//	TranscriptUtils.reAlignExtra = false;
		//	TranscriptUtils.findPolyA = false;
		//Outputs.writeGFF = false;
			IdentityProfile1.min_first_last_exon_length = 10;
			Outputs.calcBreaks=true;
			
			TranscriptUtils.coronavirus = false;
			IdentityProfile1.extra_threshold1 = 1000000;
			//Outputs.writeUnSplicedFastq = false;
			TranscriptUtils.checkAlign = false;
			TranscriptUtils.writeAnnotP = false;
			System.err.println("running in host mode");
			//CigarHash2.subclusterBasedOnStEnd = false;
			calcBreaks = false;
			SequenceUtils.mm2_splicing = "-uf";
		}
		sequential = cmdLine.getBooleanVal("sequential");
		sorted = cmdLine.getBooleanVal("sorted");
		
		String[] gffThresh_1 = cmdLine.getStringVal("gffThresh").split(":");
		Outputs.firstIsTranscriptome = cmdLine.getBooleanVal("firstIsTranscriptome");
		Outputs.library = new File(cmdLine.getStringVal("library"));
		Outputs.gffThresh = new int[gffThresh_1.length];
		
		//Outputs.gffThreshTranscriptSum=0;
		
			for(int k=0; k<gffThresh_1.length; k++){
				Outputs.gffThresh[k] = Integer.parseInt(gffThresh_1[k]);
			}
		
		
			errorAnalysis(bamFiles, RNA,reference, annotFile,readList,annotationType, 
				resDir,pattern, qual, bin, breakThresh, startThresh, endThresh,maxReads,  
				calcBreaks, filterBy5_3, annotByBreakPosition, anno, chrs, chrsToIgnore,  isoformDepthThresh, coverageDepthThresh, probInclude, fastq, chromsToRemap==null ? null: chromsToRemap.split(":"));
	}
 public static boolean sorted = true;
 static double tme0;
	public static void main(String[] args1) throws IOException, InterruptedException {
		tme0 = System.currentTimeMillis();
		CommandLine cmdLine = new ViralTranscriptAnalysisCmd2();
		String[] args = cmdLine.stdParseLine(args1);
		String bamFile = cmdLine.getStringVal("bamFile");
		String fastqFile = cmdLine.getStringVal("fastqFile");
		
		
		//Outputs.writeBed  = cmdLine.getBooleanVal("writeBed");
		Outputs.writeGFF = cmdLine.getBooleanVal("writeGFF");
		Outputs.writeIsoforms = cmdLine.getBooleanVal("writeIsoforms");

		SequenceUtils.mm2_threads = cmdLine.getIntVal("mm2_threads");
		SequenceUtils.mm2_mem = cmdLine.getStringVal("mm2_mem");
		SequenceUtils.mm2_path = cmdLine.getStringVal("mm2_path");
		SequenceUtils.mm2Preset = cmdLine.getStringVal("mm2Preset");
		String reference = cmdLine.getStringVal("reference");
		File refFile = new File(reference);
		if(!refFile.exists()) throw new RuntimeException("ref file does not exist");
		if(fastqFile!=null){
			try{
				//make a minimap index
			mm2_index = SequenceUtils.minimapIndex(refFile,  false);
			}catch(Exception exc){
				exc.printStackTrace();
			}
		}
		
		String resdir = cmdLine.getStringVal("resdir");
		String annot_file = cmdLine.getStringVal("annotation");
		
		String chroms= cmdLine.getStringVal("chroms_to_include");
		String chroms_ignore= cmdLine.getStringVal("chroms_to_ignore");
		
		//boolean gff = annot_file !=null && (annot_file.contains(".gff") || annot_file.contains(".gtf"));
		
		File resDir = new File(resdir);
		if(!resDir.exists())resDir.mkdir();
		printParams(resDir, args1);
		String[] bamFiles_ = null;
		boolean fastq= true;
		String inputFile;
		if(bamFile!=null){
			bamFiles_=bamFile.split(":");
			fastq = false;
			inputFile = bamFile;
			sorted= true;
		}else{
			
			bamFiles_=fastqFile.replace(":/", "|/").split(":");
			for(int k=0; k<bamFiles_.length; k++){
				bamFiles_[k] = bamFiles_[k].replace("|/",":/");
			}
			fastq = true;
			inputFile = fastqFile;
			sorted = false;
		}
		final boolean fastq_ = fastq;
		if(inputFile.equals("all") || inputFile.equals(".")){
			bamFiles_ = (new File("./")).list(new FilenameFilter(){
				@Override
				public boolean accept(File dir, String name) {
					return !fastq_ && name.endsWith(".bam") || (fastq_ && (name.endsWith(".fq.gz") || name.endsWith(".fastq.gz")  || name.endsWith(".fastq") || name.endsWith(".fq")));
				}
				
			});
		}
		/*String annotToRem = cmdLine.getStringVal("annotToRemoveFromGFF");
		if(annotToRem!=null){
				inner: for(int j=0; j<bamFiles_.length; j++){
					if(bamFiles_[j].indexOf(annotToRem)>=0){
						CigarCluster.annotationToRemoveFromGFF=j;
						break inner;
					}
				}
		}*/
		run(cmdLine, bamFiles_, resdir, annot_file==null ? null : new File(annot_file),  chroms, chroms_ignore, fastq, reference);
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
	
	
	/**
	 * Error analysis of a bam file. Assume it has been sorted
	 */
	static void errorAnalysis(String[] bamFiles_,boolean[] RNA, String refFile, String annot_file, String[] readList,    String annotationType, String resdir, String pattern, int qual, int round, 
			int break_thresh, int startThresh, int endThresh, int max_reads, 
			boolean calcBreaks , boolean filterBy5_3, boolean annotByBreakPosition,File gffFile, String chrToInclude, String chrToIgnore, 
			int[] writeIsoformDepthThresh, int writeCoverageDepthThresh, double probInclude, boolean fastq, String[] chromsToRemap) throws IOException {
		boolean cluster_reads = true;
		CigarHash2.round = round;
		Annotation.tolerance = round;
		
		// Integer[] basesStart = new Integer[] {0,0,0,0};
		// Integer[] basesEnd  = new Integer[] {0,0,0,0};
	//	 int leng = 10; int thresh_A = 4;int thresh_T = 5; // how many As or Ts required for determination of strand
	//	 String chars_forward = "ACGT";
	//	 String chars_reverse = "TGCA";
		 
		// Integer[] counts = new Integer[] {0,0,0};
		IdentityProfile1.annotByBreakPosition = annotByBreakPosition;
		CigarHash.cluster_by_annotation =true;// cluster_by_annotation;
		TranscriptUtils.startThresh = startThresh;
		TranscriptUtils.endThresh = endThresh;
		IdentityProfile1.writeIsoformDepthThresh =writeIsoformDepthThresh;
		IdentityProfile1.writeCoverageDepthThresh =writeCoverageDepthThresh;
		File annotSummary = new File(resdir, "annotation.csv.gz");
		if(annotSummary.exists()) annotSummary.delete();
		PrintWriter annotation_pw = new PrintWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(annotSummary, false))));
		Map<String, Integer> reads= new HashMap<String, Integer>();
		if(readList.length>0 && readList[0].length()>0){
			 Map<String, Collection<String>> map = new HashMap<String, Collection<String>>();
			 int readind =0;
			 int orfind = -1;
			 int chromind = -1;
			 int posind = -1;
			
			 for(int i=0; i<readList.length; i++){
				 InputStream is = new FileInputStream(new File(readList[i]));
				 if(readList[i].endsWith(".gz")) is = new GZIPInputStream(is);
			 BufferedReader br = new BufferedReader(new InputStreamReader(is));
				String st;
				if(readList[i].indexOf("readToCluster.txt.gz")>=0){
					st = br.readLine();
				
					List<String> head = Arrays.asList(st.split("\\t"));
					readind = head.indexOf("readID");
					orfind = head.indexOf("ORFs");
					chromind = head.indexOf("chrom");
					posind = head.indexOf("end");
				}
				while((st = br.readLine())!=null){
					String[] str  = st.split("\\s+");
					String readId = str[readind];
					String orfID = orfind>=0 && orfind < str.length ? ((chromind>=0  && chromind<str.length? str[chromind]+";" : "")+str[orfind]) : i+"";
					if(posind>=0) orfID = orfID+":"+str[posind];
					Collection<String> l= map.get(orfID) ;
					if(l==null) {
						map.put(orfID, new HashSet<String>());
						l= map.get(orfID) ;
					}
					l.add(readId);
				}
				br.close();
			 }
			readList = map.keySet().toArray(new String[0]);
			for(int i=0; i<readList.length; i++){
				Iterator<String> it =  map.get(readList[i]).iterator();
				while(it.hasNext()){
					reads.put(it.next(),i);
				}
			}
		 
		 
		}
		if(reads.size()==0){
			readList = null;
			reads=null;  //make it null to turn off this feature
		}
		IdentityProfile1.break_thresh = break_thresh;
		int len = bamFiles_.length;
		// len = 1;
		String[] in_nmes  = new String[len];
		ZipFile  anno  = null;
		boolean writeDirect = true;
		if(gffFile!=null && (gffFile.getName().indexOf(".gff")>=0 || gffFile.getName().indexOf(".gtf")>=0)){
			if(gffFile.getName().endsWith(".zip")){
				anno = new ZipFile(gffFile);
				System.err.println(anno.getName());
			}
			else {
				String out_nme = gffFile.getName();
				int ind = out_nme.lastIndexOf('.');
				out_nme = out_nme.substring(0, ind);
				File outzip = new File(gffFile.getParentFile(),out_nme+".zip");
				if(outzip.exists()){
					System.err.println("reading existing gff zip file "+outzip.getAbsolutePath());
					anno = new ZipFile(outzip);
				}
				else{
					System.err.println("making gff.zip file");
					ZipGFF gffin =  new ZipGFF( gffFile, outzip,writeDirect);
					gffin.run();
					anno = new ZipFile(outzip);
				}
			}
		}
		System.err.println("reading genomes " + (System.currentTimeMillis()- tme0)/1000+" secs");
		ArrayList<Sequence> genomes = SequenceReader.readAll(refFile, Alphabet.DNA());
		System.err.println("done reading genomes" + (System.currentTimeMillis()- tme0)/1000+" secs");
		Map<String, Integer>chrNameToIndex = new HashMap<String, Integer>();
		for(int i=0; i<genomes.size(); i++){
				chrNameToIndex.put(genomes.get(i).getName(),i);
		}
		
		// get the first chrom
		File resDir =  new File(resdir);
		if(!resDir.exists()) resDir.mkdir();
		else if(!resDir.isDirectory()) throw new RuntimeException(resDir+"should be a directory");
		///ArrayList<IdentityProfile1> profiles = new ArrayList<IdentityProfile1>();
		
		
		for (int ii = 0; ii < len; ii++) {
			String bamFile = bamFiles_[ii];
			File bam = new File( bamFile);
			in_nmes[ii] = bam.getName().split("\\.")[0];
		}
			
		
	//	genes_all_pw.close();
		
		final Iterator<SAMRecord>[] samIters = new Iterator[len];
		SamReader[] samReaders = new SamReader[len];
		boolean allNull = true;
	
		inner: for (int ii = 0; ii < len; ii++) {
			String bamFile = bamFiles_[ii];
			SamReaderFactory.setDefaultValidationStringency(ValidationStringency.SILENT);
			if(bamFile.endsWith(".bam")){
				File bam = new File( bamFile);
				samReaders[ii] = SamReaderFactory.makeDefault().open(bam);
				samIters[ii] = samReaders[ii].iterator();
			}else{
				try{
				
				samIters[ii] = SequenceUtils.getSAMIteratorFromFastq(bamFile, mm2_index, maxReads);
				}catch(Exception exc){
					System.err.println(exc.getMessage());
					System.err.println(exc.getStackTrace());
					throw new RuntimeException("could not process "+bamFile);
					//continue inner;
				}
				
			}
			allNull = false;
		}
//		if(allNull){
	//		throw new RuntimeException("No input files available "+Arrays.asList(bamFiles_));
	//	}
	//	Map<Integer, int[]> chrom_indices_to_include = new HashMap<Integer, int[]>();
		Map<String, int[]> chromsToInclude = new HashMap<String, int[]>();
		if(chrToInclude!=null && !chrToInclude.equals("all")){
			String[] chri = chrToInclude.split(":");
			for(int i=0; i<chri.length; i++){
				String[] str = chri[i].split(",");
				int[] vals ;
				if(str.length==1) vals = new int[] {0,Integer.MAX_VALUE};
				else{
					vals = new int[] {Integer.parseInt(str[1]), Integer.parseInt(str[2])};
				}
				chromsToInclude.put(str[0], vals);
			}
		}
		
		if(chromsToInclude.size()==0){
			for(int k=0; k<genomes.size(); k++){
				chromsToInclude.put(genomes.get(k).getName(), new int[] {0,Integer.MAX_VALUE});
			}
		}
		
		if(chrToIgnore!=null && chrToIgnore.equals("none")){
			String[] chri = chrToIgnore.split(":");
			for(int i=0; i<chri.length; i++){
				chromsToInclude.remove(chri[i]);
			}
		}
		
		
		/*for(int key=0; key<genomes.size();key++){
			if( chromsToInclude.containsKey(genomes.get(key).getName()))
				chrom_indices_to_include.put(key, chromsToInclude.get(key));
		}
		System.err.println("chrom indices to include");*/
		//System.err.println(chrom_indices_to_include);
		//if(chrom_indices_to_include.size()==0) chrom_indices_to_include=null;
		Collection<String> chromToRemap = Arrays.asList(chromsToRemap);
		
		
		System.err.println(chromsToInclude.size());
			int currentIndex = -1;
			
		
			
			Sequence chr = null; //genomes.get(currentIndex);
		//	Sequence chr5prime =null;Sequence chr3prime = null;
			FastqWriter[][] fqw = null;//chromToRemap.contains(chr.getName()) ? Outputs.getFqWriter(chr.getName(), resdir) : null;
				// first row is primary , second is supplementary
			
			Outputs 	outp = new Outputs(resDir,  in_nmes,   true, CigarCluster.recordDepthByPosition); 
			
			IdentityProfileHolder profile = null;
			//boolean updateProfile = true;
			
			
			Set<String> doneChr = new HashSet<String>();
			
			int numReads = 0;

			int numNotAligned = 0;
			int numSecondary =0;
			Iterator<SAMRecord> samIter= SequenceUtils.getCombined(samIters, sorted, sequential);
			float time0 = System.currentTimeMillis();//- tme0)/1000.0
			int cnt0=0;
			int ij=0;
			String previousRead = "";
			int previousRefIndex=-1;
			outer: for ( ij=0; samIter.hasNext() ;ij++ ) {
				final SAMRecord sam=samIter.next();
			    if(sam==null) break outer;
				
				if (sam.getReadUnmappedFlag()) {
					numNotAligned++;
					continue;
				}
				
				if(chromsToInclude!=null && !chromsToInclude.containsKey(sam.getReferenceName())
						&& (chromToRemap !=null && !chromToRemap.contains(sam.getReferenceName())
						)){
					System.err.println("do not have appropriate reference"+sam.getReferenceName());
					continue;
				}
				int poolID = -1;
				if(reads!=null){ 
					Integer readInd = reads.get(sam.getReadName());
					if(readInd==null) {
						if(limit_to_read_list)	continue;
					}else poolID = readInd;
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
			//	
				
				if(sam==null) break outer;
				int source_index = (Integer) sam.getAttribute(SequenceUtils.src_tag);
				//if(source_index==null) source_index = 0;
			//	System.err.println(source_index);
				if (pattern != null && (!sam.getReadName().contains(pattern)))
					continue;
				
				
				String sa  = sam.getReadString();
				boolean negFlag = sam.getReadNegativeStrandFlag();
				
			//	if(sam.getReadNegativeStrandFlag()){
			//		TranscriptUtils.flip(sam, false); // switch read but keep flag. This corrects minimaps correction when it maps
			//	}

				boolean flip = false;
				if(RNA[source_index]  ){
					flip=false;
				}
				else if(Annotation.enforceStrand){
					//at this stage cannot enforce strand with cDNA data
					throw new RuntimeException("need to implement");
					//}
				}
				
				if(flip){
					// this time we adjust flag too because we "fixing the read"  however, as yet we dont have a good way to detect if it should be flipped here
					//negFlag = !negFlag;
					
					System.err.println("flipping "+negFlag);
					TranscriptUtils.flip(sam, true);
					negFlag = sam.getReadNegativeStrandFlag();
					//sam.setReadNegativeStrandFlag(negFlag);
				}
				/*if(negFlag){
					counts[1] = counts[1] +1;
				}else{
					counts[0] = counts[0]+1;
							
				}*/
				Sequence readSeq = new Sequence(Alphabet.DNA(), sa, sam.getReadName());
				//String poolID = readSeq.get
				if (readSeq.length() <= 1) {
					//LOG.warn(sam.getReadName() +" ignored");
					
					// TODO: This might be secondary alignment, need to do something about it
					continue;
				}

				numReads++;
				if(numReads>max_reads) {
					System.err.println(numReads+" "+max_reads);
					break outer;
				}
				

				//int flag = sam.getFlags();

				if (sam.getMappingQuality() < qual) {
					numNotAligned++;
					continue;
				}
				if(probInclude<1.0 && Math.random()>probInclude){
					//randomly exclude
					continue;
				}
				
				
				
				
				// int refPos = sam.getAlignmentStart() - 1;//convert to 0-based index
				String refname = sam.getReferenceName();
				Integer refIndex = refname==null ? sam.getReferenceIndex() : chrNameToIndex.get(refname);
				if(refIndex==null) refIndex=sam.getReferenceIndex();
				// if move to another chrom, get that chrom
				if (refIndex != currentIndex) {
					System.err.println("new chrom "+refIndex + ((double)System.currentTimeMillis()- tme0)/1000.0+" secs");
					//make sure all current threads finished
					if( profile!=null && currentIndex>=0){
						{
							double timeperread = (double) ( System.currentTimeMillis() - time0)/(double) (ij-cnt0);
							System.err.println(timeperread+" milliseconds per read at "+ij);
							time0 = System.currentTimeMillis();
							cnt0= ij;
						}
						IdentityProfileHolder.waitOnThreads(100);

						profile.printBreakPoints();
						profile.getConsensus();
						doneChr.add(chr.getName());
						System.err.println("finished "+chr.getName());
					//	Outputs.waitOnThreads(100); // now we have to wait for outputs to be finished
						if(chromsToInclude != null ){
							chromsToInclude.remove(chr.getName());
							if(chromsToInclude.size()==0) {
								System.err.println("finished sample "+in_nmes[source_index]);
								break outer;
							}
						}
					}
					currentIndex = refIndex;
					
					String prev_chrom = chr==null ? "null": chr.getName();
					if(currentIndex>=genomes.size() || currentIndex <0){
						System.err.println("for some reason the refIndex greater than genome size, so finishing here "+currentIndex+" "+genomes.size());
						break outer;
					}
					chr = genomes.get(currentIndex);
					
					outp.updateChrom(chr,currentIndex );
					if(doneChr.contains(chr.getName())){
						try{
							throw new RuntimeException("not sorted contains "+ chr.getName());
						}catch(Exception exc){
							exc.printStackTrace();
						}
					}
					int seqlen = chr.length();
					Annotation annot  = null;
					if(gffFile!=null && (gffFile.getName().indexOf(".gff")>=0 || gffFile.getName().indexOf(".gtf")>=0)){
						String chrn = refname;
							annot = new GFFAnnotation(anno,chrn, seqlen, annotation_pw, gffFile.getName().indexOf(".gff")<0);
							
					}else{
						annot = annot_file == null ? new EmptyAnnotation(chr.getName(), chr.getDesc(), seqlen, annotation_pw) : 
							new Annotation(new File(annot_file), currentIndex+"", seqlen, annotation_pw, len);
					}
					boolean calcBreaks1 = calcBreaks && break_thresh < seqlen;
					
					
					 profile =  	
							new IdentityProfileHolder(genomes, chr, outp,  in_nmes, calcBreaks1, currentIndex, annot);
					
					
					
					
					if(fqw!=null){
						Outputs.waitOnThreads(Outputs.fastQwriter, 100); //make sure the fastqWriter has finished. 
						for(int i=0; i<fqw.length; i++) {
							for(int j=0; j<fqw[i].length; j++) fqw[i][j].close();
						}
					}
					fqw = chromToRemap.contains(chr.getName()) ? Outputs.getFqWriter(chr.getName(), resdir, in_nmes) : null;
				//	chr5prime = TranscriptUtils.coronavirus ? chr.subSequence(0	, Math.min( primelen, chr.length())) : null;
				//	chr3prime = TranscriptUtils.coronavirus ? chr.subSequence(Math.max(0, chr.length()-primelen), chr.length()) : null;
					System.err.println("switch chrom "+prev_chrom+"  to "+chr.getName());
				}
		
				
				if(fqw!=null){
					
						
					// this assumes that minimap2 corrected the strand and we need to reverse complement neg strand to get it back to original 
						boolean negStrand = sam.getReadNegativeStrandFlag();
						String sequence = sam.getReadString();
						int sam_ind = sam.isSecondaryOrSupplementary() ? (sam.isSecondaryAlignment() ? 1:2) : 0;
	
						/*boolean polyA = false;
						if(!negStrand && Annotation.enforceStrand){
							SWGAlignment polyAlign =  SWGAlignment.align(new Sequence(Alphabet.DNA(),sequence, sam.getReadName()), IdentityProfile1.polyA);
							if(polyAlign.getIdentity() > 0.9 * polyAlign.getLength()  && polyAlign.getLength()>15){
								polyA = true;
							}
						}
						if(polyA) sam_ind = 3;
						*/
						final FastqWriter fqw_i = fqw[sam_ind][source_index];
						Outputs.fastQwriter.execute(new Runnable(){
							public void run(){
							String baseQL = sam.getBaseQualityString();
							FastqRecord fqr= new FastqRecord(sam.getReadName(),
									negStrand ? SequenceUtil.reverseComplement(sequence): sequence,
											"",
											negStrand ?(new StringBuilder(baseQL)).reverse().toString(): baseQL
													);
							((BasicFastqWriter)fqw_i).write(fqr);
						}
					});
					continue outer;
				}
				boolean supplementary = false;
				//System.err.println(sam.getMappingQuality());
				if(sam.isSecondaryOrSupplementary()) {
					boolean supp = sam.getSupplementaryAlignmentFlag();
					if(supp && sam.getReadName().equals(previousRead) && sam.getReferenceIndex()==previousRefIndex){
						supplementary = true;
					int qual1 = sam.getMappingQuality();
					
					if(qual1 < supplementaryQ ){
						//System.err.println("remove "+qual1);
						continue; // need secondary read quality of 5
					}else{
						System.err.println("include supplementary read " + sam.getReadName()+" " + qual1);
					}
						
					}else{
					numSecondary++;
					
					continue;
					}
				}
					
				
					try{
						
						String pool = readList==null || readList.length==0 ||  poolID<0 ? null : (readList[poolID]+pool_sep);
						if(ViralTranscriptAnalysisCmd2.limit_to_read_list) pool = null;
				//		boolean flag = sam.getReadNegativeStrandFlag();
						profile.identity1(readSeq, sam, source_index, cluster_reads,  pool, q1, supplementary);
					}catch(NumberFormatException exc){
						System.err.println(readSeq.getName());
						exc.printStackTrace();
					}
					previousRead = sam.getReadName();
					previousRefIndex = sam.getReferenceIndex();
			}//samIter.hasNext()
			{
				double timeperread = (double) ( System.currentTimeMillis() - time0)/(double) (ij-cnt0);
				System.err.println(timeperread+" milliseconds per read at "+ij);
				time0 = System.currentTimeMillis();
				//cnt0= ij;
			}
			
			for(int ii=0; ii<samReaders.length; ii++){
			if(samReaders[ii] !=null) samReaders[ii].close();
			}
			
		
		
			if(profile!=null){
				IdentityProfileHolder.waitOnThreads(100);
				profile.printBreakPoints();
				profile.getConsensus();
				
			}
			if(outp!=null) outp.close();
			if(fqw!=null) {
				Outputs.waitOnThreads(Outputs.fastQwriter, 100); //make sure the fastqWriter has finished. 

				for(int i=0; i<fqw.length; i++) {
					for(int j=0; j<fqw[i].length; j++) fqw[i][j].close();
				}
			}
				
				
	if(anno!=null) anno.close();
		if(annotation_pw!=null) annotation_pw.close();
	}
}
 
