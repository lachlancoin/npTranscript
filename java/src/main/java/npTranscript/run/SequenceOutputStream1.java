package npTranscript.run;

import java.io.File;
import java.io.FileFilter;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.Stack;

import japsa.seq.Alphabet;
import japsa.seq.Sequence;
import japsa.seq.SequenceReader;
import npTranscript.cluster.Outputs;

public class SequenceOutputStream1 {
	
	
	//public static boolean keepOriginalName  = false;
	public static double se_thresh = 1.0;
	public static int tolerance = 5;
	
	public  void writeFasta(Sequence seq) throws IOException{
		
	}
	
	public static  class DetailsComparator implements Comparator<Detail>{
		 String type = "st";
		 
		 DetailsComparator(String c){
			 type = c;
		 }

		@Override
		public int compare(Detail o1, Detail o2) {
			return o1.val(type).compareTo(o2.val(type));
		}
		
		 
	 }
	
	public static void main(String[] args){
		try{
			SequenceOutputStream1.se_thresh = Double.parseDouble(args[1]);
			SequenceOutputStream1.tolerance = Integer.parseInt(args[2]);
			int minseqs = Integer.parseInt(args[3]);
			int maxseqs = Integer.parseInt(args[4]);
			SequenceOutputStream1.max_seqs_per_cluster = maxseqs;
		//	SequenceOutputStream1.keepOriginalName = true;
		//	Pattern p = Pattern.compile(args[0]);
			File[] f = new File[] {new File(args[0])};
			if(!f[0].exists()){
				 f = (new File(".")).listFiles(new FileFilter(){

					@Override
					public boolean accept(File pathname) {
						String nme = pathname.getName();
						return nme.endsWith(".fa") && nme.indexOf(args[0])>=0;
					}
					
				});
				
			}
				for(int i=0 ; i<f.length; i++){
					trim(f[i], minseqs, true);
				}
			
			
		}catch(Exception exc){
			exc.printStackTrace();
		}
	}
	 public static File trim(File target, int min_seqs, boolean filter) throws IOException {
	//	this.so = null;
		ArrayList<Sequence> genomes = SequenceReader.readAll(target.getAbsolutePath(), Alphabet.DNA());
		if(genomes.size()<min_seqs) {
			target.deleteOnExit();
			return null;
		}
		if(!filter) return target;

	/*	if(!keepOriginalName){
			target.deleteOnExit();
			target = new File(target.getParentFile(),genomes.size()+"_"+target.getName());
		}*/
    
		
		
    	
		
		
    	Detail[] det = new Detail[genomes.size()];
		
		for(int i=0; i<genomes.size(); i++){
			String[] desc = genomes.get(i).getDesc().split("\\s+");
			if(desc.length<2){
				System.err.println("WARNING, NO HEADER INFORMATION IN THE READ "+genomes.get(i).getName()+" "+target);
			}
			String[] refbreaks = desc[1].split(",");
			det[i] = new Detail(i);
			int start = Integer.parseInt(refbreaks[0]);
			int end= Integer.parseInt(refbreaks[refbreaks.length-1]);
			det[i].st= start;
			det[i].en = end;
			if(refbreaks.length>2){
				int br_start = Integer.parseInt(refbreaks[1]);
				int br_end = Integer.parseInt(refbreaks[refbreaks.length-2]);
				det[i].break_L = br_start;// - start;
				det[i].break_R = br_end;
			}else{
				//det[i].span_l = 0;//end;//-start;
			}
		}
		/*; 
		Arrays.sort(en);
		*/
		
		double[][] meanSE = Detail.calcMeanSD(det);
		
		int mid = (int) Math.floor((double) genomes.size()*perc_target);
		Arrays.sort(det, new DetailsComparator("st"));
		int st_target = det[mid].st; 
		Arrays.sort(det, new DetailsComparator("en"));
		int en_target = det[det.length-mid].en;
		Arrays.sort(det, new DetailsComparator("index"));

		SequenceOutputStream1 so = new SequenceOutputStream1(target);
		int excl_count =0;
		for(int i=0; i<genomes.size(); i++){
			if(det[i].index!=i) throw new RuntimeException("!!");
	
			Sequence seq = genomes.get(i);
			int start = det[i].st;
			int end= det[i].en;
			boolean include = det[i].include(meanSE[0], meanSE[1],se_thresh, tolerance);
			if(start<st_target+10 && end>en_target - 10 && include){
			int read_st_new = Math.max(0, st_target-start); //truncate if st[i] < st_target 
			int read_end_new = seq.length()-Math.max(0,  end-en_target ); //truncate if en[i]> en_target 
		//	System.err.println(read_st_new+" "+read_end_new+" "+seq.length());
			
			if(read_st_new < read_end_new){
				Sequence seq1 = seq.subSequence(read_st_new, read_end_new);
				seq1.setName(seq.getName());
				seq1.setDesc(seq.getDesc()+" "+read_st_new+","+read_end_new+","+seq1.length());
				so.stack.push(seq1);
			}else{
				excl_count++;
				//System.err.println("excluded "+seq.getName()+" from clusters");
			}
			}else{
				excl_count++;
//				System.err.println("excluded "+seq.getName()+" from clusters");
			}
		}
		
		if(excl_count>0) System.err.println("excluded "+excl_count+" of "+genomes.size());
		if(so.stack.size()>min_seqs){
			so.printAll();
		}else{
			System.err.println(target.getName()+" did not meet minimum record count after trimming " + so.stack.size()+" < "+min_seqs);
			so.stack.clear();
		}
		//so.close();
		return so.target;
	}
	
	static double perc_target = 0.5; // more means greater truncation
  public static int max_seqs_per_cluster = 100000;
 private  File target ;
  //boolean append;
  int thresh = 10;
  int seqs_printed=0;
  
  public SequenceOutputStream1(File out) {
	
	this.target = out;
	//this.append = append;
	}

  

  /*public SequenceOutputStream1(ZipOutputStream outS) {
	so = new OutputStreamWriter(outS);
  }*/


   boolean lock = false;  // is stream locked for closing
public void printAll() throws IOException {
	lock = true; 
	Runnable run = new Runnable(){
		@Override
		public void run() {
		//	System.err.println("launching "+target);
			try{
			 OutputStreamWriter	so = new OutputStreamWriter(new FileOutputStream(target, true));
 			while(stack.size()>0 && seqs_printed < max_seqs_per_cluster){
				
				 Sequence seq = stack.pop();
				 String[] sequ1 = new String[] {seq.toString()};
					String name = seq.getName();
					String desc = seq.getDesc();
					boolean allnull = true;
					for(int i=0; i<sequ1.length; i++){
						if(sequ1[i]!=null) allnull=false;
					}
					if(allnull) continue;
					so.write(">"+name+" "+desc+"\n");
					int step = Integer.MAX_VALUE;
					for(int k=0; k<sequ1.length; k++){
						if(sequ1[k]!=null){
						for(int j=0; j<sequ1[k].length(); j+=step){
							String substr = sequ1[k].substring(j, Math.min(sequ1[k].length(), j+step));
							//if(writeNumb) os1.write(j+" ");
							so.write(substr);
							so.write("\n");
						}
						}
					}
					seqs_printed++;
				 
				 
				
				
			}
 			stack.clear();
		
			 so.close();
			//	System.err.println("done.. "+target+" ");
			}catch(IOException exc){
				exc.printStackTrace();
			}
			lock = false;
			
		}
		
	};
	Outputs.writeCompressDirsExecutor.execute(run);
}

  
  public void write(Sequence seq) throws IOException {
	  if(seqs_printed < max_seqs_per_cluster){
		stack.push(seq);
		if(stack.size()==thresh){
			this.printAll();
		}
	  }
	}

Stack<Sequence>  stack= new Stack();


/*public void flush() throws IOException{
	if(so!=null) {
		if(lock) close = true; else so.close();
	}
	// TODO Auto-generated method stub
	
}*/


/*
public boolean close()  throws IOException{
	if(so!=null) {
		if(lock) close = true; else so.close();
		return true;
	}else{
		return false;
	}
	// TODO Auto-generated method stub
	
}*/



public String entryname() {
	// TODO Auto-generated method stub
	return target.getName();
}






  
 
  
  
}
