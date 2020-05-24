package npTranscript.run;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.util.Stack;
import java.util.zip.ZipOutputStream;

import japsa.seq.Sequence;
import japsa.seq.SequenceOutputStream;

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
