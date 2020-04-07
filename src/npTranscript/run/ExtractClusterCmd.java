package npTranscript.run;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

import japsa.util.CommandLine;
import japsa.util.deploy.Deployable;

@Deployable(scriptName = "npTranscript.run", scriptDesc = "Clustering of results from ViralTranscriptAnalysisCmd2")
public class ExtractClusterCmd extends CommandLine {
	
	
	public ExtractClusterCmd() {
		super();
		Deployable annotation = getClass().getAnnotation(Deployable.class);
		setUsage(annotation.scriptName() + " [options]");
		setDesc(annotation.scriptDesc());

		addString("inDir", null, "Name of inputDir", true);
		addStdHelp();
	}
	
	
	
	
	public static void main(String[] args) throws IOException, InterruptedException {
		CommandLine cmdLine = new ExtractClusterCmd();
		args = cmdLine.stdParseLine(args);

		String inDir = cmdLine.getStringVal("inDir");
		File in = new File(inDir);
		File transcriptF = new File(in, "0transcripts.txt.gz");
		File outdir = new File(in, "clusters");
		outdir.mkdir();
		File readsF = new File(in, "0readToCluster.txt.gz");
	
		
		try{
			ExtractClusterCmd ec = new ExtractClusterCmd();
			ec.run(transcriptF, readsF, outdir);
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
	Table transcripts;
	File outdir;
	
	
	
	public void run(File transcriptsF , File readsF,  File outdir) throws IOException{
		reads = new Table(readsF,"\t");
		transcripts = new Table(transcriptsF,"\t");
	    this.outdir = outdir;
		transcripts.sort("countTotal", true);
		transcripts.filter("type_nme",  "5_3");
		transcripts.limit("countTotal", 10, true);
		String[] ids = transcripts.getCol("ID");
		Map<String, String> m= new HashMap<String, String>();
		for(int i=0; i<ids.length; i++){
			List<String> all_ids = Arrays.asList(ids[i].split(";"));
			Table t2 = this.reads.extract("clusterId", "upstream", "downstream", all_ids);
			t2.print(outdir);
			for(Iterator<String> it = all_ids.iterator(); it.hasNext();){
				m.put(it.next(), t2.name);
			}
			t2.print(outdir);
		}
		transcripts.replaceAll("ID", m);
		transcripts.print(outdir);
		
	}
	

		
		
	}


