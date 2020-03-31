package npTranscript.cluster;

import org.apache.commons.math3.Field;
import org.apache.commons.math3.FieldElement;
import org.apache.commons.math3.exception.MathArithmeticException;
import org.apache.commons.math3.exception.NullArgumentException;


/**
 * @author Lachlan Coin
 *
 */

public class IntegerField implements FieldElement<Integer>{
	Integer x;
	IntegerField(int x) {
		this.x = x;
	}
	public Integer getVal(){
		return x;
	}
	@Override
	public Integer add(Integer a) throws NullArgumentException {
		
		return  x+a;
	}

	@Override
	public Integer subtract(Integer a) throws NullArgumentException {
		return x - a;
	}

	@Override
	public Integer negate() {
		return -x;
	}

	@Override
	public Integer multiply(int n) {
		return x*n;
	}

	@Override
	public Integer multiply(Integer a) throws NullArgumentException {
		return x*a;
	}

	@Override
	public Integer divide(Integer a) throws NullArgumentException, MathArithmeticException {
		return null;
	}

	@Override
	public Integer reciprocal() throws MathArithmeticException {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public Field<Integer> getField() {
		// TODO Auto-generated method stub
		return null;
	}

	
}