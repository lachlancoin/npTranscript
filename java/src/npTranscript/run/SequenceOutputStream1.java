package npTranscript.run;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Stack;
import java.util.TreeSet;
import java.util.zip.ZipOutputStream;

import japsa.seq.Alphabet;
import japsa.seq.Sequence;
import japsa.seq.SequenceReader;

public class SequenceOutputStream1 {
	
	public  void writeFasta(Sequence seq) throws IOException{
		String[] sequ1 = new String[] {seq.toString()};
		String name = seq.getName();
		String desc = seq.getDesc();
		boolean allnull = true;
		for(int i=0; i<sequ1.length; i++){
			if(sequ1[i]!=null) allnull=false;
		}
		if(allnull) return;
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
	}
	
	public void trim(int min_seqs) throws IOException {
		this.so = null;
		ArrayList<Sequence> genomes = SequenceReader.readAll(target.getAbsolutePath(), Alphabet.DNA());
		target.delete();
    
		
		if(genomes.size()<min_seqs) return;
    	target = new File(target.getParentFile(),genomes.size()+"_"+target.getName());

		int[] st = new int[genomes.size()];
		int[] en = new int[genomes.size()];
		for(int i=0; i<genomes.size(); i++){
			String[] desc = genomes.get(i).getDesc().split("\\s+");
			String[] refbreaks = desc[1].split(",");
			
			int start = Integer.parseInt(refbreaks[0]);
			int end= Integer.parseInt(refbreaks[refbreaks.length-1]);
			st[i] = start;
			en[i] = end;
		}
		Arrays.sort(st); Arrays.sort(en);
		
		int mid = (int) Math.floor((double) genomes.size()*perc_target);
		int st_target = st[mid]; int en_target = en[en.length-mid];
		for(int i=0; i<genomes.size(); i++){
			String[] desc = genomes.get(i).getDesc().split("\\s+");
			String[] refbreaks = desc[1].split(",");
		//	String[] readbreaks = desc[2].split(",");
			Sequence seq = genomes.get(i);
			int start = Integer.parseInt(refbreaks[0]);
			int end= Integer.parseInt(refbreaks[refbreaks.length-1]);
			int read_st_new = Math.max(0, st_target-start); //truncate if st[i] < st_target 
			int read_end_new = seq.length()-Math.max(0,  end-en_target ); //truncate if en[i]> en_target 
		//	System.err.println(read_st_new+" "+read_end_new+" "+seq.length());
			if(read_st_new < read_end_new){
				Sequence seq1 = seq.subSequence(read_st_new, read_end_new);
				seq1.setName(seq.getName());
				seq1.setDesc(seq.getDesc()+" "+read_st_new+","+read_end_new+","+seq1.length());
				this.stack.push(seq1);
			}else{
				System.err.println("excluded "+seq.getName()+" from clusters");
			}
		}
		if(this.stack.size()>min_seqs){
			this.printAll();
		}else{
			this.stack.clear();
		}
		this.close();
		
	}
	static double perc_target = 0.5; // more means greater truncation
  private OutputStreamWriter so; 
  File target ;
  boolean append;
  int thresh = 10;
  
  public SequenceOutputStream1(File out, boolean append) {
	
	this.target = out;
	this.append = append;
	}

  

  public SequenceOutputStream1(ZipOutputStream outS) {
	so = new OutputStreamWriter(outS);
  }


public void printAll() throws IOException {
	if(so==null) so = new OutputStreamWriter(new FileOutputStream(target, append));
	while(stack.size()>0){
		writeFasta(stack.pop());;
	}
}

  
  public void write(Sequence seq) throws IOException {
	  if(so!=null) writeFasta(seq);
	  else{
		stack.push(seq);
		if(stack.size()==thresh){
			this.printAll();
		}
	  }
	}

Stack<Sequence>  stack= new Stack();


public void flush() throws IOException{
	if(so!=null) so.close();
	// TODO Auto-generated method stub
	
}



public boolean close()  throws IOException{
	if(so!=null) {
		so.close();
		return true;
	}else{
		return false;
	}
	// TODO Auto-generated method stub
	
}



public String entryname() {
	// TODO Auto-generated method stub
	return target.getName();
}






  
 
  
  
}
