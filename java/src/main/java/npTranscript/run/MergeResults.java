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
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.io.PrintWriter;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.SortedMap;
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;

import ch.systemsx.cisd.hdf5.HDF5Factory;
import ch.systemsx.cisd.hdf5.IHDF5Reader;
import ch.systemsx.cisd.hdf5.IHDF5Writer;
import japsa.util.CommandLine;
import japsa.util.deploy.Deployable;

/**
 * @author Lachlan Coin
 *
 */

@Deployable(scriptName = "npTranscript.run", scriptDesc = "Analysis of coronavirus sequence data")
public class MergeResults extends CommandLine {

	public MergeResults() {
		super();
		Deployable annotation = getClass().getAnnotation(Deployable.class);
		setUsage(annotation.scriptName() + " [options]");
		setDesc(annotation.scriptDesc());
		addString("dirs", null, "Input h5", false);
		addString("match", null, "Input h5", false);
		addInt("thresh",10,"threshold",false);
		addString("toInclude",null, "which transcripts to include", false);
	}
	 
	public static Map<String, String> decodeM = new HashMap<String, String>();
	public static void main(String[] args1) throws IOException, InterruptedException {
		CommandLine cmdLine = new MergeResults();
		/*File f = new File("./");
		
		/File[] dirs = f.listFiles(new FileFilter(){

			@Override
			public boolean accept(File pathname) {
				return pathname.isDirectory() ;
			}
			
		});*/
		String[] args = cmdLine.stdParseLine(args1);
	//	String h5file = cmdLine.getStringVal("h5File");
		String[] dirs1 = cmdLine.getStringVal("dirs").split(":");
		File[] dirs = new File[dirs1.length];
		String[] match = cmdLine.getStringVal("match").split(":");
		thresh = cmdLine.getIntVal("thresh");
		for(int i=0; i<dirs.length; i++){
			dirs[i] = new File(dirs1[i]);
		//	inds[i] 
		}
		
		File decode = new File("decode.txt");
		if(decode.exists()){
			BufferedReader br = new BufferedReader(new FileReader(decode));
			String st = "";
			while((st = br.readLine())!=null){
				String[] str = st.split("\\s+");
				decodeM.put(str[0], str[1]);
			}
			br.close();

		}
		
		File outdir = new File("combined");
		outdir.mkdir();
		String toIncl = cmdLine.getStringVal("toInclude");
		List<String> toInclude =toIncl==null ? null : Arrays.asList(toIncl.split(":"));
		//Arrays.fill(inds, 0);
		Set<String> to_exclude= new HashSet<String>();
		Set<String> excluded = new HashSet<String>();
		int[] inds = mergeIsoforms("0.isoforms.h5", dirs, false, "/trans:/counts".split(":"), match, outdir, to_exclude, excluded, toInclude);
		to_exclude.addAll(excluded);
		int[] inds1 = mergeIsoforms("0.clusters.h5", dirs, true, "/depth:/depthStart:/depthEnd".split(":"), match, outdir, to_exclude, excluded, toInclude);
		System.err.println("size "+to_exclude.size());
		mergeText(dirs, "0.annot.txt.gz",  inds, outdir);
		//String[] tocop  
		//Files.copy(source, target, StandardCopyOption.REPLACE_EXISTING);
		//Files.copy(, target, options)

		

	}
	public static void mergeText(File[] dirs, String in, int[] inds, File out) throws FileNotFoundException, IOException{
		BufferedReader[] br  = new BufferedReader[dirs.length];
		String[][] head = new String[br.length][];
		StringBuffer sb1 = new StringBuffer();
		for(int i=0; i<dirs.length; i++){
			br[i] = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(new File(dirs[i], in)))));
			head[i] = br[i].readLine().split("\t");
		
		}
		sb1.append(head[0][0]+"\t"+head[0][1]+"\t"+head[0][2]);
		for(int i=0; i<dirs.length;i++){
			sb1.append("\t"+"Spliced_5_"+i+"\tUnspliced_5_"+i);
		}
		PrintWriter pw = new PrintWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(new File(out,in)))));
		pw.println(sb1.toString());
		int start = 3;
		
		
		outer: while(true){
			StringBuffer sb = new StringBuffer();
			
			for(int i=0; i<br.length; i++){
				int st1 = inds[i]*2+start;
				String str = br[i].readLine();
				if(str==null) break outer;
				List<String> l = Arrays.asList(str.split("\t"));
				if(i==0){
					sb.append(l.get(0)+"\t"+l.get(1)+"\t"+l.get(2));
				}
				sb.append("\t"+l.get(st1)+"\t"+l.get(st1+1)); //need 
			}
			pw.println(sb.toString());
		}
		for(int i=0; i<br.length; i++) br[i].close();
		pw.close();
		
	}
	public static int thresh = 10;
	public static void merge(IHDF5Reader[] reader,IHDF5Writer writer,String loc, int[] inds, boolean array, int offset,
			Set<String> to_exclude, Set<String> excluded, List<String> toInclude){ //loc=/trans
		int[] cnts = new int[reader.length];
		HashSet<String>ins = new HashSet<String>();
		HashSet<String>[]ins_i = new HashSet[reader.length];
		for(int i=0; i<reader.length; i++){
			ins_i[i] = new HashSet<String>(reader[i].getGroupMembers(loc));
			ins.addAll(ins_i[i]);
		}
		SortedMap[] obj = new SortedMap[reader.length];
		for(Iterator<String> it = ins.iterator(); it.hasNext();){
			String nxt = it.next();
			String location=loc+"/"+nxt;
			if(to_exclude.contains(nxt)) {
				System.err.println("skpping "+location);
				continue;
			}
			if(toInclude!=null && !toInclude.contains(nxt)) {
				continue;
			}
			Arrays.fill(cnts,0);
			if(array){
				for(int i=0; i<reader.length; i++){
					if(ins_i[i].contains(nxt)){
						int[][] cnts_i = reader[i].readIntMatrix(location);
						obj[i] = convertToHash(cnts_i, inds[i], offset);
					}else{
						obj[i] = new TreeMap<Integer, int[]>();
					}
				}
				int[][] merged = merge(obj, inds);//, cnts);
				//int sumcnt=sum(cnts);
				//if(sumcnt>thresh)				{
					writer.writeIntMatrix(location, merged);
			//	}
			}else{
				
				for(int i=0; i<reader.length; i++){
					if(ins_i[i].contains(nxt)){
						int[] cnts_i = reader[i].readIntArray(location);
						cnts[i] = cnts_i[inds[i]];
					}
				}
				int sumcnt=sum(cnts);
				if(sumcnt>thresh){
				writer.writeIntArray(location, cnts);
				}else{
					System.err.println("excluding "+nxt);
					excluded.add(nxt);
				}
			}
		//	System.err.println("h");
		}
	}
	private static int sum(int[] cnts) {
		int sum=0;
		for(int i=0; i<cnts.length; i++) sum+=cnts[i];
		return sum;
	}


	private static SortedMap<Integer, int[]> convertToHash(int[][] cnts_i, int ind, int offset) {
		SortedMap<Integer, int[] >m = new TreeMap<Integer, int[]>();
		for(int i=0; i<cnts_i.length; i++){
			int[] v = cnts_i[i];
			int len = v.length-offset;
			if(len%2 !=0) throw new RuntimeException("!!");
			int len1 = (len)/2;
		
			int[] v1 = new int[] {v[2*ind+offset],v[2*ind+offset+1]}; // needs work for multi index files
			m.put(v[0], v1);
		}
		return m;
	}


	private static int[][] merge(SortedMap<Integer, int[]>[] obj, int[] inds) {
		// TODO Auto-generated method stub
		TreeSet<Integer> s = new TreeSet<Integer>();
		for(int i=0; i< obj.length ;i++){
			s.addAll(obj[i].keySet());
		}
		int len1 = obj.length;
		int start1 = 1;
		int start2 = 1+obj.length;
		int len = obj.length*2 + 1;
		int[][] res = new int[s.size()][len];//k,
		int k=0;
		for(Iterator<Integer> it = s.iterator(); it.hasNext();k++){
			Integer key = it.next();
			int[] v = res[k];
			Arrays.fill(v,0);
			v[0] = key;
			for(int i=0; i<obj.length; i++){
				if(obj[i].containsKey(key)){
					int[] v1 = obj[i].get(key);
					v[2*i+1] = v1[0];
					v[2*i+2] = v1[1];
					//cnts[i]  = Math.max(v1[0], cnts[i]);
				}
			}
			
		}
		return res;
	}


	public static int[] mergeIsoforms(String h5file, File[] dirs, boolean array, String[] address, String[] match, File out,
			Set<String> toExclude, Set<String> excluded, List<String> toInclude
			){
		File[] h5 = new File[dirs.length];
		int offset = 0;
		if(array) offset=1;  // to account for the pos
		IHDF5Reader[] reader=  new IHDF5Reader[h5.length];
		String[]	sources =  new String[reader.length+offset];
		if(offset==1) sources[0] = "pos";
		int[] counts  = new int[reader.length];
		int[] inds = new int[dirs.length];
		
		for(int i=0; i<reader.length; i++){
			File f= new File(dirs[i],h5file);
			
			reader[i] = HDF5Factory.openForReading(f.getAbsolutePath());
			String[] header = reader[i].readStringArray("header");
			inds[i] = Arrays.asList(header).indexOf(match[i])-offset;
			if(inds[i]<0)inds[i]=0;
			sources[i+offset] = header[inds[i]+offset];
			String decodeSt = decodeM.get(sources[i+offset]);
			if(decodeSt!=null){
				sources[i+offset]=decodeSt;
			}
		}
		
		File outf = new File(out,h5file);
		outf.delete();
		IHDF5Writer writer = HDF5Factory.open(outf);
		writer.writeStringArray("header", sources);
		for(int k=0; k<address.length; k++){
			merge(reader, writer,address[k], inds, array, offset, toExclude, excluded, toInclude);
		}
		System.err.println("done");
		writer.close();
		return inds;
	}

}
