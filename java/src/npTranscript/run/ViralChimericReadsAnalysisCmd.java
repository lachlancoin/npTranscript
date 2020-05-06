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
		addInt("hitThresh", 200, "Length threshold for hits on the genome");
		addInt("leftOverThresh", 200, "Length threshold for left over sequences");
		addStdHelp();
	}
	
	public static void main(String[] args1) {
		CommandLine cmdLine = new ViralChimericReadsAnalysisCmd();
		String[] args = cmdLine.stdParseLine(args1);
		String bamFile = cmdLine.getStringVal("bamFile");
		String outPrefix = cmdLine.getStringVal("outPrefix");
		hit_thresh = cmdLine.getIntVal("hitThresh");
		leftover_thresh = cmdLine.getIntVal("leftOverThresh");
		run(bamFile.split(":"), outPrefix);
	}
	
	private static int hit_thresh = 200;
	private static int leftover_thresh = 200;
	
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
	
	public static void run(String[] files, String out_prefix) {
		String out_fq = out_prefix+"_complete.fastq";
		String out_leftover = out_prefix+"_leftover.fastq";
		String out_bam = out_prefix+".bam";
		
		if(new File(out_fq).exists()||new File(out_leftover).exists()||new File(out_bam).exists()) {
			throw new RuntimeException("Output file exsits!!!");
		}
		
		try {
			BufferedWriter bw_fq = new BufferedWriter(new FileWriter(new File(out_fq)));
			BufferedWriter bw_leftover = new BufferedWriter(new FileWriter(new File(out_leftover)));
			
			SAMFileWriter samWriter = null;
			for(String f : files) {
				List<String> bam_files = new ArrayList<>();
				File file = new File(f);
				if(!file.exists()) {
					bw_fq.close();
					bw_leftover.close();
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

						StringBuilder descBuilder = new StringBuilder();
						for(Block block : blocks) {
							descBuilder.append(block.q1);
							descBuilder.append(",");
							descBuilder.append(block.q2);
							descBuilder.append(",");
							descBuilder.append(block.s1);
							descBuilder.append(",");
							descBuilder.append(block.s2);
							descBuilder.append(",");
							descBuilder.append(block.strand);
							descBuilder.append(";");
						}
						descBuilder.setLength(descBuilder.length()-1);
						String desc = descBuilder.toString();
						
						RangeSet<Integer> leftover = TreeRangeSet.create();
						leftover.add(Range.closed(1, readLen).canonical(DiscreteDomain.integers()));
						//ImmutableRangeSet<Integer> leftover = ImmutableRangeSet.copyOf(fullRange);
						leftover.removeAll(rangeCov);
						
						boolean chimera = false;
						int k = 1;
						for(Range<Integer> r : leftover.asRanges()) {
							int st = r.lowerEndpoint();
							int ed = r.upperEndpoint()-1;
							String seqstr = sequence.substring(st, ed);
							if(seqstr.length()>=leftover_thresh&&!isPolyATailorRepeat(seqstr)) {
								String qualstr = qual.substring(st, ed);
								bw_leftover.write("@"+read_id+"."+k+" "+st+","+ed+","+readLen+" "+desc+"\n");
								bw_leftover.write(seqstr+"\n");
								bw_leftover.write("+\n");
								bw_leftover.write(qualstr+"\n");
								chimera = true;
								++k;
							}
						}
						if(chimera) {
														
							bw_fq.write("@"+read_id+" "+desc+"\n");
							bw_fq.write(sequence+"\n");
							bw_fq.write("+\n");
							bw_fq.write(qual+"\n");
							
							for(SAMRecord r : records) samWriter.addAlignment(r);	
						}

						if(record!=null) read_id = record.getReadName();
					}

					iterator.close();
					samReader.close();
				}
			}
			
			samWriter.close();
			bw_fq.close();
			bw_leftover.close();
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
		return max+secmax>(b[0]+b[1]+b[2]+b[3]+b[4])*0.8;
	}
}












