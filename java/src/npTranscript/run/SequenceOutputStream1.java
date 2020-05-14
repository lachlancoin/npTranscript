package npTranscript.run;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.util.Stack;
import java.util.zip.ZipOutputStream;

import japsa.seq.Sequence;
import japsa.seq.SequenceOutputStream;

public class SequenceOutputStream1 {
  private SequenceOutputStream so; 
  File target ;
  boolean append;
  int thresh = 10;
  
  public SequenceOutputStream1(File out, boolean append) {
	
	this.target = out;
	this.append = append;
	}

  

  public SequenceOutputStream1(ZipOutputStream outS) {
	so = new SequenceOutputStream(outS);
  }


public void printAll() throws IOException {
	if(so==null) so = new SequenceOutputStream(new FileOutputStream(target, append));
	while(stack.size()>0){
		stack.pop().writeFasta(so);;
	}
}

  
  public void write(Sequence seq) throws IOException {
	  if(so!=null) seq.writeFasta(so);
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
