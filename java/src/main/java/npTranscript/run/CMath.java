package npTranscript.run;

import com.sun.jna.Library;
import com.sun.jna.Native;
import com.sun.jna.Platform;

public interface CMath extends Library {
	  double cosh(double value);
	  CMath INSTANCE = Native.load(Platform.isWindows() ? "msvcrt" : "c", CMath.class);
	  
	  public static void main(String[] args){
		  //CMath lib = Native.load(Platform.isWindows()?"msvcrt":"c", CMath.class);
		  
		  double result = INSTANCE.cosh(0);
		  System.err.println(result);
	  }
}
