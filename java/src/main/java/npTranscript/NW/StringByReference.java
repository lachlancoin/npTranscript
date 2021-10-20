package npTranscript.NW;

import com.sun.jna.ptr.ByReference;

public class StringByReference extends ByReference {

/*	protected StringByReference(int dataSize) {
		super(dataSize);
		// TODO Auto-generated constructor stub
	}*/
	 public StringByReference() {
	        super(4);
	    }
	 
	 public StringByReference(String value) {
	        super(40);
	        setValue(value);
	    }

	    public void setValue(String value) {
	        getPointer().setString(0, value);
	    }

	    public String getValue() {
	        return getPointer().getString(0);
	    }
		int length() {
			return getValue().length();
		}

	   /* @Override
	    public String toString() {
	        return String.format("int@0x%1$x=0x%2$x (%2$d)", Pointer.nativeValue(getPointer()), getValue());
	    }*/

}
