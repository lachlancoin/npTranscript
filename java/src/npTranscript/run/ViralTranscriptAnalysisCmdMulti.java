package npTranscript.run;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Date;
import java.util.Map;

import japsa.seq.JapsaAnnotation;
import japsa.util.CommandLine;
import npTranscript.cluster.GFFAnnotation;

public class ViralTranscriptAnalysisCmdMulti extends ViralTranscriptAnalysisCmd2 {
  
	
	
	public static void main(String[] args1) throws IOException, InterruptedException {
		CommandLine cmdLine = new ViralTranscriptAnalysisCmd2();
		String[] args = cmdLine.stdParseLine(args1);
		String[] bamFile = cmdLine.getStringVal("bamFile").split(":");
		String[] resdir = cmdLine.getStringVal("resdir").split(":");
		File[] resDir = new File[resdir.length];
		for(int i=0; i<resdir.length; i++){
			resDir[i] =  new File(resdir[i]);
			if(!resDir[i].exists()) resDir[i].mkdir();
			printParams(resDir[i], args1);
		}
		String annot_file = cmdLine.getStringVal("annotation");
		boolean gff = annot_file.contains(".gff");
		Map<String, JapsaAnnotation> anno  = null;
		boolean SARS = !gff;
		if(gff){
			String  annotationType = cmdLine.getStringVal("type");
			System.err.println("reading annotation");
			 anno =  GFFAnnotation.readAnno(annot_file, annotationType);
			System.err.println("done");
		}
		for(int i=0; i<bamFile.length; i++){
			run(cmdLine, bamFile[i], resdir[i], anno,SARS);
		}
		
		
		
		
	
	}
}
