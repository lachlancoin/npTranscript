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
 * 05/05/2020 - Chenxi Zhou                                       
 ****************************************************************************/

package npTranscript.run;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.FilenameFilter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;

import com.google.common.collect.DiscreteDomain;
import com.google.common.collect.Range;
import com.google.common.collect.RangeSet;
import com.google.common.collect.TreeRangeSet;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.DefaultSAMRecordFactory;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.util.SequenceUtil;
import japsa.util.CommandLine;
import japsa.util.deploy.Deployable;

/**
 * @author Chenxi Zhou
 *
 */

@Deployable(scriptName = "npChimera.run", scriptDesc = "Detection of chimeric reads")
public class ViralChimericReadsAnalysisCmd extends CommandLine {
	
	public ViralChimericReadsAnalysisCmd() {
		super();
		Deployable annotation = getClass().getAnnotation(Deployable.class);
		setUsage(annotation.scriptName() + " [options]");
		setDesc(annotation.scriptDesc());
		addString("bamFile", null, "Name of bam file", true);
		addString("outPrefix", null, "Prefix for output files", true);
		addString("reference", null, "Reference file in FASTA format", true);
		addInt("hitThresh", 200, "Length threshold for hits on the genome");
		addInt("clipThresh", 20, "Lenght threshold for clipping.");
		addInt("leftOverThresh", 200, "Length threshold for left over sequences");
		addInt("flankSize", 20, "Flanking size for the reference sequence at fusion point.");
		addStdHelp();
	}
	
	public static void main(String[] args) {
		CommandLine cmdLine = new ViralChimericReadsAnalysisCmd();
		cmdLine.stdParseLine(args);
		String bamFile = cmdLine.getStringVal("bamFile");
		String outPrefix = cmdLine.getStringVal("outPrefix");
		String reference = cmdLine.getStringVal("reference");
		hit_thresh = cmdLine.getIntVal("hitThresh");
		clip_thresh = cmdLine.getIntVal("clipThresh");
		leftover_thresh = cmdLine.getIntVal("leftOverThresh");
		if(clip_thresh>leftover_thresh) {
			System.err.println("WRANING: leftOverThresh is smaller than clipThresh. "
					+ "Set leftOverThresh="+clip_thresh);
			leftover_thresh = clip_thresh;
		}
		flank_size = cmdLine.getIntVal("flankSize");
		run(bamFile.split(":"), reference, outPrefix);
	}
	
	private static int hit_thresh;
	private static int leftover_thresh;
	private static int clip_thresh;
	private static int flank_size;
	
	private static class Block {
		private String str_id;
		private int s1, s2, q1, q2;
		private String strand;

		public Block(String str_id, int s1, int s2, int q1, int q2, String strand) {
			this.str_id = str_id;
			this.s1 = s1;
			this.s2 = s2;
			this.q1 = q1;
			this.q2 = q2;
			this.strand = strand;
		}
	}
	
	public static void run(String[] files, String reference, String out_prefix) {
		String out_fq = out_prefix+"_complete.fastq";
		String out_leftover = out_prefix+"_leftover.fastq";
		String out_bam = out_prefix+".bam";
		
		if(new File(out_fq).exists()||new File(out_leftover).exists()||new File(out_bam).exists()) {
			throw new RuntimeException("Output file exsits!!!");
		}
		
		try {
			BufferedWriter bw_complete = new BufferedWriter(new FileWriter(new File(out_fq)));
			BufferedWriter bw_leftover = new BufferedWriter(new FileWriter(new File(out_leftover)));
			IndexedFastaSequenceFile dict = new IndexedFastaSequenceFile(new File(reference));
				
			SAMFileWriter samWriter = null;
			for(String f : files) {
				List<String> bam_files = new ArrayList<>();
				File file = new File(f);
				if(!file.exists()) {
					bw_complete.close();
					bw_leftover.close();
					dict.close();
					if(samWriter!=null) samWriter.close();
					throw new RuntimeException("!!!");
				}
				if(file.isFile()) {
					bam_files.add(f);
				} else if(file.isDirectory()) {
					File dir = new File(f);
					File[] fs = dir.listFiles(new FilenameFilter() {
						public boolean accept(File dir, String name) {
							return name.toLowerCase().endsWith(".bam");
						}
					});
					for(File s : fs)
						bam_files.add(s.getAbsolutePath());
				}
				
				for(String bam_file : bam_files) {
					SamReader samReader = SamReaderFactory.make()
							.validationStringency(ValidationStringency.LENIENT)
							.samRecordFactory(DefaultSAMRecordFactory.getInstance())
							.open(new File(bam_file));
					if(samWriter==null)
						samWriter = new SAMFileWriterFactory().makeSAMOrBAMWriter(samReader.getFileHeader(),
				                true, new File(out_bam));
					SAMRecordIterator iterator = samReader.iterator();
					List<SAMRecord> records = new ArrayList<>();
					List<Block> blocks = new ArrayList<>();
					SAMRecord record = iterator.hasNext()?iterator.next():null;
					String read_id = record==null?null:record.getReadName();
					
					while( record!=null ) {
						if(record.getReadUnmappedFlag()) {
							record = iterator.hasNext()?iterator.next():null;
							read_id = record==null?null:record.getReadName();
							continue;
						}
						records.clear();
						if(!record.isSecondaryAlignment()&&!isPolyATailorRepeat(record)) 
							records.add(record);
						while((record = iterator.hasNext()?iterator.next():null)!=null &&
								record.getReadName().equals(read_id)) {
							if(!record.isSecondaryAlignment()&&!isPolyATailorRepeat(record))
								records.add(record);
						}
						
						RangeSet<Integer> rangeCov = TreeRangeSet.create();
						blocks.clear();
						
						int readLen = 0;
						String sequence = null;
						String qual = null;
						for(SAMRecord r : records) {
							if(!r.isSecondaryOrSupplementary()) {
								readLen = r.getReadLength();
								sequence = r.getReadString();
								qual = r.getBaseQualityString();
								if(r.getReadNegativeStrandFlag()) {
									sequence = SequenceUtil.reverseComplement(sequence);
									qual = new StringBuilder(qual).reverse().toString();
								}
							}

							int hclip = 0;
							CigarElement c = r.getCigar().getFirstCigarElement();
							if(c.getOperator()==CigarOperator.HARD_CLIP)
								hclip = c.getLength();

							int s1 = r.getAlignmentStart();
							int s2 = r.getAlignmentEnd();
							int q1 = r.getReadPositionAtReferencePosition(s1)+hclip;
							int q2 = r.getReadPositionAtReferencePosition(s2)+hclip;

							if(r.getReadNegativeStrandFlag()) {
								int tmp = q1;
								q1 = readLen-q2+1;
								q2 = readLen-tmp+1;
							}
							
							rangeCov.add(Range.closed(q1, q2).canonical(DiscreteDomain.integers()));
							
							blocks.add(new Block(r.getReferenceName(), s1, s2, q1, q2, r.getReadNegativeStrandFlag()?"r":"f"));
						}
						
						int alnLen = 0;
						for(Range<Integer> r : rangeCov.asRanges()) {
							alnLen += r.upperEndpoint()-r.lowerEndpoint();
						}
						
						if(readLen==0||alnLen<hit_thresh) { 
							//primary alignment was filtered out or alignment too short
							if(record!=null) read_id = record.getReadName();
							continue;
						}
						
						Collections.sort(blocks, new Comparator<Block>() {

							@Override
							public int compare(Block o1, Block o2) {
								// TODO Auto-generated method stub
								return Integer.compare(o1.q1, o2.q1);
							}

						});
						
						// 1. check if the alignment chunks are consistent (for supplementary alignments)
						//     - alignment strand need to be consistent for all alignment chunks
						//     - no gaps or overlap over 20bp on read sequence between multiple alignment chunks
						//     - no overlap over 20bp on reference between multiple alignment chunks (gaps are allowed)
						// 2. fusion in only one side
						//     - if the clip at 5'- or 3'-end is over 20bp, will check if it's polyA tail or not
						//       if not will consider as a clip
						//     - clip is only allowed in one side
						// 3. output of leftover sequences
						//     - if the clip sequence is over 200bp, will output as leftover sequences
						//       leftover sequence will be only at 5'-end or 3'-end (otherwise will be filtered out by #2)
						//     - the reference sequence flanking 20bp of the upstream and downstream at the fusion point
						//       will be always 41bp by default - padded with 'N' if flanking sequence is too short
						// 4. the FASTQ header of the leftover sequences will contain three information fields
						//     - add '-y' option for alignment with minimap2
						//     - SAMRecord getStringAttribute() could be used to retrieve these information
						//       record.getStringAttribute("AP") for the alignment position of the original sequence
						//       record.getStringAttribute("LO") for the position of the leftover sequence
						//       record.getStringAttribute("FS") for the flanking sequence on the reference sequence
						
						Block block = blocks.get(0);
						String refId = block.str_id;
						boolean isPositiveStrand = block.strand.equals("f");
						int qstart = block.q1;
						int qend   = block.q2;
						int sstart = isPositiveStrand?block.s1:block.s2;
						int send   = isPositiveStrand?block.s2:block.s1;
						boolean badAlignment = false;
						
						for(int i=1; i<blocks.size(); i++) {
							block = blocks.get(i);
							int olap = isPositiveStrand ? (send-block.s1) : (block.s2-send);
							if(!block.str_id.equals(refId) ||
									isPositiveStrand!=block.strand.equals("f") ||
									Math.abs(qend-block.q1)>20 ||
									olap > 20) {
								badAlignment = true;
								break;
							}
							qend = block.q2;
							send = isPositiveStrand ? block.s2 : block.s1;
						}
						
						boolean leftOverAt5End = qstart>clip_thresh&&!isPolyATail(sequence.substring(0, qstart-1));
						boolean leftOverAt3End = qend<=readLen-clip_thresh&&!isPolyATail(sequence.substring(qend));
						
						if(badAlignment||leftOverAt5End==leftOverAt3End) {
							if(record!=null) read_id = record.getReadName();
							continue;
						}

						// build describe line
						StringBuilder descBuilder = new StringBuilder("AP:Z:");
						for(Block b : blocks) {
							descBuilder.append(b.q1);
							descBuilder.append(",");
							descBuilder.append(b.q2);
							descBuilder.append(",");
							descBuilder.append(b.s1);
							descBuilder.append(",");
							descBuilder.append(b.s2);
							descBuilder.append(",");
							descBuilder.append(b.strand);
							descBuilder.append(";");
						}
						descBuilder.setLength(descBuilder.length()-1);
						String desc = descBuilder.toString();
						
						String refSequence = dict.getSequence(refId).getBaseString();
						boolean isChimera = false;
						if(qstart>leftover_thresh){
							String seqstr = sequence.substring(0, qstart-1);
							String qualstr = qual.substring(0, qstart-1);
							String refstr = getFlankSequence(refSequence, sstart-1, flank_size);
							if(!isPositiveStrand) refstr = SequenceUtil.reverseComplement(refstr);
							bw_leftover.write("@"+read_id+" LO:Z:"+1+","+(qstart-1)+","+readLen+" "+desc+" "+"FS:Z:"+refstr+"\n");
							bw_leftover.write(seqstr+"\n");
							bw_leftover.write("+\n");
							bw_leftover.write(qualstr+"\n");
							isChimera = true;
						}
						if(qend<=readLen-leftover_thresh) {
							String seqstr = sequence.substring(qend);
							String qualstr = qual.substring(qend);
							String refstr = getFlankSequence(refSequence, send-1, flank_size);
							if(!isPositiveStrand) refstr = SequenceUtil.reverseComplement(refstr);
							bw_leftover.write("@"+read_id+" LO:Z:"+(qend+1)+","+readLen+","+readLen+" "+desc+" "+"FS:Z:"+refstr+"\n");
							bw_leftover.write(seqstr+"\n");
							bw_leftover.write("+\n");
							bw_leftover.write(qualstr+"\n");
							isChimera = true;
						}
						
						if(isChimera) {		
							bw_complete.write("@"+read_id+" "+desc+"\n");
							bw_complete.write(sequence+"\n");
							bw_complete.write("+\n");
							bw_complete.write(qual+"\n");
							
							for(SAMRecord r : records) samWriter.addAlignment(r);	
						}

						if(record!=null) read_id = record.getReadName();
					}

					iterator.close();
					samReader.close();
				}
			}
			
			samWriter.close();
			bw_complete.close();
			bw_leftover.close();
			dict.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	private static boolean isPolyATailorRepeat(SAMRecord record) {
		// TODO Auto-generated method stub
		String seqstr = record.getReadString();
		if(!record.isSecondaryOrSupplementary()) {
			int readLen = record.getReadLength();
			int st = 0, ed = 0;
			CigarElement c = record.getCigar().getFirstCigarElement();
			if(c.getOperator().isClipping()) st = c.getLength();
			c = record.getCigar().getLastCigarElement();
			if(c.getOperator().isClipping()) ed = c.getLength();
			ed = readLen-ed;
			seqstr = seqstr.substring(st, ed);
		}
		
		return isPolyATailorRepeat(seqstr);
	}
		
	private static boolean isPolyATailorRepeat(String seqstr) {	
		int[] b = new int[5];
		for(char c : seqstr.toCharArray()) {
			switch(Character.toUpperCase(c)) {
				case 'A':
					b[0]++;
					break;
				case 'T':
					b[1]++;
					break;
				case 'U':
					b[2]++;
					break;
				case 'C':
					b[3]++;
					break;
				case 'G':
					b[4]++;
					break;
				default:
					break;	
			}
		}
		
		int max=0, secmax=0;
		for(int i : b) {
			if(i>=max) {
				secmax = max;
				max = i;
			} else if(i>secmax) {
				secmax = i;
			}
		}
		return max+secmax>=(b[0]+b[1]+b[2]+b[3]+b[4])*0.8;
	}
	
	private static boolean isPolyATail(String seqstr) {	
		int[] b = new int[3];
		for(char c : seqstr.toCharArray()) {
			switch(Character.toUpperCase(c)) {
				case 'A':
					b[0]++;
					break;
				case 'T':
					b[1]++;
					break;
				case 'U':
					b[2]++;
					break;
				default:
					break;	
			}
		}
		
		double max = seqstr.length()*0.8;
		return b[0]>=max||b[1]>=max||b[2]>=max;
	}
	
	private static String getFlankSequence(String sequence, int p, int flank) {
		// TODO Auto-generated method stub
		StringBuilder flank_str = new StringBuilder();
		int readLen = sequence.length();
		flank_str.setLength(0);
		if(p<flank) {
			for(int k=0; k<flank-p; k++)
				flank_str.append('N');
		}
		flank_str.append(sequence.substring(Math.max(0, p-flank), Math.min(readLen, p+flank+1)));
		if(p+flank+1>readLen) {
			for(int k=0; k<p+flank+1-readLen; k++)
				flank_str.append('N');
		}
		return flank_str.toString();
	}
}












