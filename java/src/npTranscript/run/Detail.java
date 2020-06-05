package npTranscript.run;

public class Detail{
	public int st,en,break_L, break_R, index;
	Detail(int index){
		this.index = index;
	}
	
	private void add(double[] mean) {
		// TODO Auto-generated method stub
		mean[0] = mean[0] + st;
		mean[1] = mean[1] + en;
		mean[2] = mean[2] + break_L;
		mean[3] = mean[3] + break_R;
	}
	private void addVar(double[] var, double[] mean) {
		// TODO Auto-generated method stub
		var[0] = var[0] + Math.pow(st-mean[0],2);
		var[1] = var[1] + Math.pow(en-mean[1],2);
		var[2] = var[2] + Math.pow(break_L-mean[2],2);
		var[3] = var[3] + Math.pow(break_R-mean[3],2);
	}
	public String toString(){
		return st+"_"+en+"_"+break_L+"_"+break_R+"_"+index;
	}
	public Integer val(String type){
		if(type.equals("st")){
			return st;
		}else if(type.equals("en")){
			return en;
		}else if(type.equals("index")){
		return index;
	}
		return null;
	}
	
	public static double[][] calcMeanSD(Detail[] det) {
		double[] mean = new double[4];
		double[] se = new double[4];
		for(int i=0; i<det.length; i++){
			det[i].add(mean);
		}
		for(int i=0; i<mean.length; i++){
			mean[i] = mean[i]/(double)det.length;
		}
		for(int i=0; i<det.length; i++){
			det[i].addVar(se,mean);
		}
		for(int i=0; i<mean.length; i++){
			se[i] = Math.sqrt(se[i]/(double)det.length);
		}
		return new double[][] {mean,se};
	}

	public boolean include(double[] ds, double[] ds2, double d, int tol) {
		double se2 = Math.max(ds2[2]*d, tol);
		double se3 = Math.max(ds2[3]*d, tol);
		return this.break_L>=ds[2]-se2 && break_L<=ds[2]+se2 && 
				this.break_R>=ds[3]-se3 && break_L<=ds[3]+se3;
	}
	
	
}