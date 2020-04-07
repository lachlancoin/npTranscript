package npTranscript.run;

import java.io.File;
import java.io.FileFilter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

import japsa.util.CommandLine;
import japsa.util.deploy.Deployable;
import npTranscript.cluster.CigarHash;

@Deployable(scriptName = "npTranscript.run", scriptDesc = "Merging results from ViralTranscriptAnalysisCmd2")
public class MergeCmd extends CommandLine {
	
	
	public MergeCmd() {
		super();
		Deployable annotation = getClass().getAnnotation(Deployable.class);
		setUsage(annotation.scriptName() + " [options]");
		setDesc(annotation.scriptDesc());

		addString("inDir","./", "Name of inputDir", false);
		addStdHelp();
	}
	
	//ID	start	end	type_nme	breaks	hash	startBreak	endBreak	leftGene	rightGene	totLen	countTotal	count0count1	depth0	depth1	errors0	errors1	error_ratio0	error_ratio1
//  	9	8	29858	5_3	0,1,2624,2989	-1287377339	69	26237	LEADER	E	0	3	3	0	11012	0	707	00.0642	NaN

	
	
	public static void main(String[] args) throws IOException, InterruptedException {
		CommandLine cmdLine = new MergeCmd();
		args = cmdLine.stdParseLine(args);
		String tf_name = "0transcripts.txt.gz";

		File inDir_ = new File(cmdLine.getStringVal("inDir"));
		File[] inDir = inDir_.listFiles(new FileFilter(){

			@Override
			public boolean accept(File pathname) {
				return pathname.isDirectory() && (new File(pathname, tf_name).exists());
			}
			
		});
		File[] transcriptF = new File[inDir.length];
		for(int i=0; i<inDir.length; i++){
			transcriptF[i] = new File(inDir[i], tf_name);
		}
	
		
		try{
			MergeCmd ec = new MergeCmd();
			ec.run(transcriptF, tf_name);
		}catch(Exception exc){
			exc.printStackTrace();
		}
	}
	
  static int cntr = 1; 
  static String split = "\t";
	static class StringComp implements Comparator<String[]>{
		final int col;
		final int mult;
		
		StringComp(int col, boolean decr){
			this.col = col;
			mult = decr ? -1 : 1;
		}
		@Override
		public int compare(String[] o1, String[] o2) {
			Integer i1 = Integer.parseInt(o1[col]);
			Integer i2 = Integer.parseInt(o2[col]);
			return mult*i1.compareTo(i2);
		}

		
		
	}
	Table reads;
	Table[] transcripts;
	File outdir = new File("merged");
	
	Map<String, String[][]> map = new HashMap<String, String[][]>();
	
	
	static class MergeInds{
		List<Integer> remaining = new ArrayList<Integer>();
		Integer[] merge_inds;
		private String[] target;
		int[] merge_offset;
		String[][] head_input;
		final String default_val;
		public String[] merge( String[][]input , boolean add_j ){
			
			boolean first= true;
			for(int j=0; j<input.length; j++){
				String[] str = input[j];
				if(str!=null && first){
					first = false;
					for(int i=0; i<remaining.size(); i++){
						target[i] = str[remaining.get(i)];
					}
				}
				for(int k=0; k<merge_inds.length; k++){
						String nme = str==null ? default_val : str[merge_inds[k]];
						if(add_j) nme = nme.replaceAll("0" , ""+j);
						target[merge_offset[k]+j] = nme;
				}
			}
			return target.clone();
		}
		public String[] newHeader(){
			return this.merge(head_input, true);
		}
		
		public MergeInds(List<String> header, String[] merge, String[] toignore, int num_sources, String default_val){
			merge_inds = new Integer[merge.length];
			merge_offset = new int[merge.length];
			this.default_val = default_val;
			for(int i=0; i<header.size(); i++){
				remaining.add(i);
			}
				
			for(int i=0; i<merge_inds.length; i++){
				merge_inds[i] = header.indexOf(merge[i]);
				if(merge_inds[i]<0) {
					throw new RuntimeException("could not find "+merge[i]);
				}
				remaining.remove(merge_inds[i]);
			}
			for(int i=0; i<toignore.length; i++){
				Integer torem = header.indexOf(toignore[i]);
				if(torem<0) throw new RuntimeException("!!");
				remaining.remove(torem);
			}
			target = new String[remaining.size() + num_sources*merge.length];
			{
			int i= remaining.size();
			for(int k=0; k<merge_inds.length; k++){
					merge_offset[k] = i;
					i+=num_sources;
			}
			}
			 head_input  = new String[num_sources][];
			for(int i=0; i<num_sources; i++){
				head_input[i] = header.toArray(new String[0]);
			}
		
		//	merge(merged_header, input, true);
			
		}
	}
	
	
	public void run(File[] transcriptsF , String tf_name) throws IOException{
		String[] ids = new String[] {"type_nme", "leftGene", "rightGene"};
		List<String> header = null;
		int num_sources = transcriptsF.length;
		StringBuffer info = new StringBuffer(); 
		for(int i=0; i<transcriptsF.length; i++){
			info.append(transcriptsF[i].getParentFile().getName());
			 info.append("\t");
			Table transcripts = new Table(transcriptsF[i],"\t");
			header = transcripts.header;
			transcripts.addToMap(map, ids, i, num_sources);
		}
		String[] merge = new String[] {"count0",  "depth0", "errors0", "error_ratio0"};
		String[] toignore = new String[] {"countTotal"};
		MergeInds mi = new MergeInds(header,  merge,  toignore, num_sources, "0");
		List<String> head_new = Arrays.asList(mi.newHeader());
		Table t1 = new Table(head_new, tf_name);
		for(Iterator<String[][] > it = map.values().iterator(); it.hasNext();){
			String[]str = mi.merge(it.next(), false);
			t1.data.add(str);
		}
		
		t1.print(outdir, info.toString());
		
	}
	

		
		
	}


