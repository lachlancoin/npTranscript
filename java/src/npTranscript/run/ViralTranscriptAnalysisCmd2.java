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
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Date;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import japsa.seq.Alphabet;
import japsa.seq.JapsaAnnotation;
import japsa.seq.Sequence;
import japsa.seq.SequenceReader;
import japsa.util.CommandLine;
import japsa.util.deploy.Deployable;
import npTranscript.cluster.Annotation;
import npTranscript.cluster.CigarHash;
import npTranscript.cluster.CigarHash2;
import npTranscript.cluster.GFFAnnotation;
import npTranscript.cluster.IdentityProfile1;
import npTranscript.cluster.IdentityProfile1.Outputs;
import npTranscript.cluster.TranscriptUtils;

/**
 * @author Lachlan Coin
 *
 */

@Deployable(scriptName = "npTranscript.run", scriptDesc = "Analysis of coronavirus sequence data")
public class ViralTranscriptAnalysisCmd2 extends CommandLine {
//	private static final Logger LOG = LoggerFactory.getLogger(HTSErrorAnalysisCmd.class);

	
	public ViralTranscriptAnalysisCmd2() {
		super();
		Deployable annotation = getClass().getAnnotation(Deployable.class);
		setUsage(annotation.scriptName() + " [options]");
		setDesc(annotation.scriptDesc());

		addString("bamFile", null, "Name of bam file", true);
		//addString("breaks", null, "Position File, for looking for specific breaks");
		addString("reference", null, "Name of reference genome", true);
		addString("annotation", null, "ORF annotation file or GFF file", true);
		addString("readList", null, "List of reads", false);
	//	addString("genesToInclude", null, "Names of genes to include in analysis, only used if annotation is GFF file. "
	//			+ " Example: --genesToInclude Name=ACE2:Name=TMPRSS2", false);
		String all_types = "gene:ncRNA_gene:pseudogene";			
			
			
				
				
			
				
				

		addString("type", all_types, "Type of annotation (only included if annotation is GFF file", false);
		addString("chroms", "all", "Restrict to these chroms, colon delimited", false);

//		String genesToInclude = "Name=ACE2:Name=TMPRSS2";
		addString("resdir", "results"+System.currentTimeMillis(), "results directory");

		addInt("maxReads", Integer.MAX_VALUE, "ORF annotation file");
		

		addString("pattern", null, "Pattern of read name, used for filtering");
		addInt("qual", 0, "Minimum quality required");
		addInt("bin", 1, "Bin size for numerical hashing");
		addInt("breakThresh", 10, "Thresh for break points to match clusters");

		addInt("startThresh", 100, "Threshold for having 5'");
		addInt("endThresh", 100, "Threshold for having 3'");
		//addDouble("overlapThresh", 0.95, "Threshold for overlapping clusters");
	//	addBoolean("coexpression", false, "whether to calc coexperssion matrices (large memory for small bin size)");
		addBoolean("overwrite", false, "whether to overwrite existing results files");
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
	
 public static void run(CommandLine cmdLine, String bamFile, String resDir,Map<String, JapsaAnnotation> anno, boolean SARS, Set<String> chrs) throws IOException{
		String reference = cmdLine.getStringVal("reference");
		int qual = cmdLine.getIntVal("qual");
		int bin = cmdLine.getIntVal("bin");
		int breakThresh = cmdLine.getIntVal("breakThresh");
		String pattern = cmdLine.getStringVal("pattern");
		String annotFile = cmdLine.getStringVal("annotation");
		String readList = cmdLine.getStringVal("readList");
		//String genesToInclude = cmdLine.getStringVal("genesToInclude");
		String  annotationType = cmdLine.getStringVal("type");
		boolean overwrite  = cmdLine.getBooleanVal("overwrite");
		int startThresh = cmdLine.getIntVal("startThresh");
		int endThresh = cmdLine.getIntVal("endThresh");
		int maxReads = cmdLine.getIntVal("maxReads");
		
		boolean sorted = true;
		boolean calcBreaks = false;// whether to calculate data for the break point heatmap, true for SARS_COV2
		boolean filterBy5_3 = false;// should be true for SARS_COV2
		boolean annotByBreakPosition = false;  // should be true for SARS_COV2
		if(SARS){
			calcBreaks  = true; 
			filterBy5_3 = true;
			annotByBreakPosition = true;
		}
			errorAnalysis(bamFile, reference, annotFile,readList,annotationType, 
				resDir,pattern, qual, bin, breakThresh, startThresh, endThresh,maxReads,  sorted , 
				calcBreaks, filterBy5_3, annotByBreakPosition, anno, chrs, overwrite);
	}
	public static void main(String[] args1) throws IOException, InterruptedException {
		CommandLine cmdLine = new ViralTranscriptAnalysisCmd2();
		String[] args = cmdLine.stdParseLine(args1);
		String bamFile = cmdLine.getStringVal("bamFile");
		String resdir = cmdLine.getStringVal("resdir");
		String annot_file = cmdLine.getStringVal("annotation");
		String chroms= cmdLine.getStringVal("chroms");
		Set<String> chrs = null;
		if(!chroms.equals("all")){
			chrs = new HashSet<String>(Arrays.asList(chroms.split(":")));
		}
		boolean gff = annot_file.contains(".gff");
		Map<String, JapsaAnnotation> anno  = null;
		boolean SARS = !gff;
		if(gff){
			String  annotationType = cmdLine.getStringVal("type");
			System.err.println("reading annotation");
			 anno =  GFFAnnotation.readAnno(annot_file, annotationType, chrs);
			// GFFAnnotation.read
			System.err.println("done reading annotation");
		}
		File resDir = new File(resdir);
		if(!resDir.exists())resDir.mkdir();
		printParams(resDir, args1);
		run(cmdLine, bamFile, resdir, anno, SARS, chrs);
		
		
		// paramEst(bamFile, reference, qual);
	}

	
	/**
	 * Error analysis of a bam file. Assume it has been sorted
	 */
	static void errorAnalysis(String bamFiles, String refFile, String annot_file, String readList,    String annotationType, String resdir, String pattern, int qual, int round, 
			int break_thresh, int startThresh, int endThresh, int max_reads,  boolean sorted,
			boolean calcBreaks , boolean filterBy5_3, boolean annotByBreakPosition,Map<String, JapsaAnnotation> anno, Set<String>chrToInclude, boolean overwrite ) throws IOException {
		boolean cluster_reads = true;
		CigarHash2.round = round;
		IdentityProfile1.annotByBreakPosition = annotByBreakPosition;
		CigarHash.cluster_by_annotation =true;// cluster_by_annotation;
		TranscriptUtils.startThresh = startThresh;
		TranscriptUtils.endThresh = endThresh;
		boolean calcTree = false;
		String[] bamFiles_ = bamFiles.split(":");
		Set<String> reads = null;
		if(readList!=null){
			reads = new HashSet<String>();
			BufferedReader br = new BufferedReader(new FileReader(new File(readList)));
			String st;
			while((st = br.readLine())!=null){
				String st_ = st.split("\\s+")[0];
				System.err.println(st_);
				reads.add(st_);
			}
			br.close();
		}
		
		
		TranscriptUtils.break_thresh = break_thresh;
		int len = bamFiles_.length;
		// len = 1;
		String[] in_nmes  = new String[len];
		ArrayList<Sequence> genomes = SequenceReader.readAll(refFile, Alphabet.DNA());

		// get the first chrom
		File resDir =  new File(resdir);
		if(!resDir.exists()) resDir.mkdir();
		else if(!resDir.isDirectory()) throw new RuntimeException(resDir+"should be a directory");
		///ArrayList<IdentityProfile1> profiles = new ArrayList<IdentityProfile1>();
		
		boolean gff = anno!=null;
		for (int ii = 0; ii < len; ii++) {
			String bamFile = bamFiles_[ii];
			File bam = new File( bamFile);
			in_nmes[ii] = bam.getName().split("\\.")[0];
		}
		Outputs outp = new Outputs(resDir,  in_nmes, overwrite); 
		Set<String> doneChr = new HashSet<String>();
			
		
	//	genes_all_pw.close();
		IdentityProfile1 profile = null;
		outer1: for (int ii = 0; ii < len; ii++) {
			int source_index = ii;
			int currentIndex = 0;
			boolean skipChrom = false;
			Sequence chr = genomes.get(currentIndex);
			//profile=  profiles.get(currentIndex);
			String bamFile = bamFiles_[ii];
			File bam = new File( bamFile);
		//	in_nmes[ii] = bam.getName().split("\\.")[0];
		//	for (int jj = 0; jj < genomes.size(); jj++) {
			//	if(profiles.get(jj)!=null){
				//	profiles.get(jj).updateSourceIndex(ii);
				//}
			//}
			
			SamReaderFactory.setDefaultValidationStringency(ValidationStringency.SILENT);
			SamReader samReader = null;// SamReaderFactory.makeDefault().open(new File(bamFile));

			if ("-".equals(bamFile))
				samReader = SamReaderFactory.makeDefault().open(SamInputResource.of(System.in));
			else
				samReader = SamReaderFactory.makeDefault().open(bam);

			SAMRecordIterator samIter = samReader.iterator();
			// Read the reference genome
			

			long totReadBase = 0, totRefBase = 0;
			int numReads = 0;

			int numNotAligned = 0;
	
			outer: for (int cntr = 0; samIter.hasNext() && cntr < max_reads; cntr++) {
				SAMRecord sam = samIter.next();
				if(readList!=null && !reads.contains(sam.getReadName())){
					cntr=cntr-1;
					continue outer;
					
				}
				if (pattern != null && (!sam.getReadName().contains(pattern)))
					continue;

				// make the read seq
				Sequence readSeq = new Sequence(Alphabet.DNA(), sam.getReadString(), sam.getReadName());
				if (readSeq.length() <= 1) {
					// LOG.warn(sam.getReadName() +" ignored");
					// TODO: This might be secondary alignment, need to do something about it
					continue;
				}

				numReads++;

				if (sam.getReadUnmappedFlag()) {
					numNotAligned++;
					continue;
				}

				int flag = sam.getFlags();

				if (sam.getMappingQuality() < qual) {
					numNotAligned++;
					continue;
				}

				// int refPos = sam.getAlignmentStart() - 1;//convert to 0-based index
				int refIndex = sam.getReferenceIndex();
				
				// if move to another chrom, get that chrom
				if (refIndex != currentIndex) {
					if( profile!=null && currentIndex>=0){
						profile.printBreakPoints();
						profile.getConsensus();
						profile= null;
						doneChr.add(chr.getName());
						System.err.println("finished "+chr.getName()+" "+in_nmes[source_index]);
						chrToInclude.remove(chr.getName());
						if(chrToInclude != null && chrToInclude.size()==0) {
							System.err.println("finished sample "+in_nmes[source_index]);
							break outer;
						}
					}
					currentIndex = refIndex;
					chr = genomes.get(currentIndex);
					
					
				}
				if(chrToInclude!=null && ! chrToInclude.contains(chr.getName())){
					skipChrom=true;
				}else{
					skipChrom = false;
				}
				if(!skipChrom){
					if(profile==null){
							if(doneChr.contains(chr.getName())){
								try{
									throw new RuntimeException("not sorted contains "+ chr.getName());
								}catch(Exception exc){
									exc.printStackTrace();
								}
							//	System.err.println("warning done"+chr.getName());
							}
							int seqlen = chr.length();
						//	if(!anno.containsKey(chr))
							Annotation annot = null;
							if(gff){
								JapsaAnnotation annot1 = anno.get(chr.getName());
								if(annot1==null){
									try{
										throw new RuntimeException("no annotation for  "+ chr.getName());
									}catch(Exception exc){
										exc.printStackTrace();
									}
								}else{
									annot = new GFFAnnotation(annot1, seqlen);
								}
							}else{
								annot = new Annotation(new File(annot_file), currentIndex+"", seqlen);
							}
							profile = new IdentityProfile1(chr, outp,  in_nmes, startThresh, endThresh, annot, calcBreaks, chr.getName(), currentIndex);
					}
				
				 
			
					try{
					TranscriptUtils.identity1(chr, readSeq, sam, profile, source_index, cluster_reads, chr.length());
					}catch(NumberFormatException exc){
						System.err.println(readSeq.getName());
						exc.printStackTrace();
					}
				}

			}
			samReader.close();
			
		}
		
			if(profile!=null){
				profile.printBreakPoints();
				profile.getConsensus();
				
			}
		
		outp.close();
		ExtractClusterCmd.cluster(resdir,outp.transcripts_file, outp.reads_file, filterBy5_3);
	
	}

}

/*
 * RST* ----------------------------------------------------------
 * jsa.hts.errorAnalysis*: Error analysis of sequencing data
 * ----------------------------------------------------------
 * 
 * jsa.hts.errorAnalysis* assesses the error profile of sequencing data by
 * getting the numbers of errors (mismatches, indels etc) from a bam file.
 * Obviously, it does not distinguish sequencing errors from mutations, and
 * hence consider mutations as errors. It is best to use with the bam file from
 * aligning sequencing reads to a reliable assembly of the sample.
 * 
 * <usage>
 * 
 * RST
 */