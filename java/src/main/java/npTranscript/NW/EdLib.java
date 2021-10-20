package npTranscript.NW;

import edlib.EdlibAlignConfig;
import edlib.EdlibAlignResult;
import edlib.EdlibLibrary;


//import edlib.*;

public class EdLib  {
	
	
	  
	  public static void main(String[] args){
		//   config =  
		   EdlibAlignConfig.ByValue config =EdlibLibrary.INSTANCE.edlibDefaultAlignConfig();
		  StringByReference query = new StringByReference("hello");
		   StringByReference target =new StringByReference( "world!");
	EdlibAlignResult res=	
		EdlibLibrary.INSTANCE.edlibAlign(query.getPointer(), 5, target.getPointer(), 6, config);
	//	INSTANCE.edlibAlign(query, 5, target, 6, config);

	System.err.println(res.editDistance);
	  }
	 
	  
	 
}
