package npTranscript.cluster;

import org.apache.commons.math3.Field;


/**
 * @author Lachlan Coin
 *
 */

public class IField implements Field<IntegerField>{
	
		IntegerField zero = new IntegerField(0);
		IntegerField one = new IntegerField(1);
		@Override
		public IntegerField getZero() {
			// TODO Auto-generated method stub
			return zero;
		}

		@Override
		public IntegerField getOne() {
			// TODO Auto-generated method stub
			return one;
		}

		@Override
		public Class getRuntimeClass() {
			// TODO Auto-generated method stub
			return IntegerField.class;
		}

		
		
	
}