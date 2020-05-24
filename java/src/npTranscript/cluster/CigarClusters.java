package npTranscript.cluster;

import java.io.IOException;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;

import japsa.seq.Sequence;


/**
 * @author Lachlan Coin
 *
 */


public class CigarClusters {
	
	
	
//	final double thresh;
	
	public CigarClusters(int num_sources){
	 this.num_sources = num_sources;
	 this.seqlen =0;
	 this.refseq =null;
	 this.annot = null;
	}
	
	private static String getKey(String[] str, int[] inds) {
		StringBuffer sb = new StringBuffer();
		for(int i=0; i<inds.length; i++) {
			sb.append(str[inds[i]]);
			sb.append(",");
		}
		return sb.toString();
	}

	CigarClusters( Sequence refSeq, Annotation annot, int num_sources){
	//	this.thresh = thresh;
		this.refseq = refSeq;
		this.seqlen = refSeq.length();
		this.annot = annot;
		this.num_sources = num_sources;
	}
	/*public DistanceMatrix getDistanceMatrix( PrintWriter pw){
		CigarCluster[] l1 = this.l.values().toArray(new CigarCluster[0]);
		int len  = l1.length;
		double[][] res = new double[len][];
		String[] labels = new String[len];
		for(int i=0; i<len; i++) {
			
			res[i] = new double[len];
			CigarCluster cc = l1[i];
			labels[i] = cc.id;
			res[i][i] =0;
			for(int j=0; j<i; j++) {
				CigarCluster cc_j = l1[j];
				double dist = 1-cc.similarity(cc_j);
				res[i][j] = dist;
				res[j][i] = dist;
			}
			
		}
		for(int i=0; i<len; i++) {
			pw.print(labels[i]+","+l1[i].index+",");
			pw.println(getString(res[i]));
		}
	
		IdGroup group = new  SimpleIdGroup(labels);
		DistanceMatrix dm = new DistanceMatrix(res, group);
		return dm;
	}*/
	Map<CigarHash, CigarCluster> l = new HashMap<CigarHash, CigarCluster>();


	
	public void matchCluster(CigarCluster c1,  int source_index, int num_sources, int chrom_index,String[] clusterIDs) throws NumberFormatException{
		String clusterID;
		CigarHash2 subID ;
		if(!l.containsKey(c1.breaks_hash)){
			CigarCluster newc = new CigarCluster("ID"+chrom_index+"."+l.keySet().size(), num_sources, c1, source_index);
			clusterID = newc.id();
			l.put(newc.breaks_hash, newc);
			subID = newc.breaks;
		}else{
			CigarCluster clust = l.get(c1.breaks_hash);
			subID = clust.merge(c1, num_sources, source_index);
			clusterID = clust.id();
		}
	clusterIDs[0] = clusterID;
	clusterIDs[1] = "_"+subID.toString(CigarHash2.round)+"_";
		
		//return clusterID;
	}

	/*private void writeSeq(SequenceOutputStream seqFasta,Annotation annot,  PrintWriter exonP, int[][] exons, CigarCluster cc, String refseq, int[] firstlast){
		StringBuffer descline; String annotline;
		String id = cc.id+"";
		int read_count = cc.readCountSum;
		if(cc.breakSt>0){
			String leftseq = refseq.subSequence(cc.start, cc.breakSt).toString();
			String rightseq = refseq.subSequence(cc.breakEnd, cc.end).toString();
		}else{
			String rightseq = refseq.subSequence(cc.start, cc.end).toString();
		}
		for(int j=0; j<exons.length; j++) {
			int start = exons[j][0];
			int end = exons[j][1];
			exonP.println(id+"\t"+start+"\t"+end+"\t"+read_count);

			
			annotline.append(annot.calcORFOverlap(start, end, first_last, transcript_len));

			int len = end-start+1;
			descline.append(";");
			descline.append(start); descline.append("-"); descline.append(end); descline.append(","); descline.append(len);
			subseq.append();
			
			transcript_len += len;
			//seqline.append(subseq.toString());
		}
		Sequence subseq1 = new Sequence(refseq.alphabet(),subseq.toString().toCharArray(), id);
	
		descline.append(" "); descline.append(annotline);
		subseq1.setDesc(descline.toString());
		subseq1.writeFasta(seqFasta); 
	}*/
	//static int tolerance = 10;
	final Annotation annot;
	final Sequence  refseq;
	final int seqlen;
	final int num_sources;
	public void process(CigarCluster cc,   Outputs o,String chrom,int chrom_index){
		cc.addZeros(seqlen); 
		String id = cc.id();
		int totalDepth = cc.readCountSum();
		if(Outputs.doMSA!=null && Outputs.doMSA.contains(cc.getTypeNme(seqlen))){
			try{
			if(o.mergeSourceClusters){
				o.msa(id,0,cc);	
			}else{
			for(int k=0; k<this.num_sources; k++){
				o.msa(id,k,cc);
				
			}
			}
				
			}catch(Exception exc){
				exc.printStackTrace();
			}
		}
		if(o.clusterW!=null && totalDepth>IdentityProfile1.writeCoverageDepthThresh){
		
		 int[][] matr =cc.getClusterDepth(num_sources, this.refseq);
		
		 o.writeIntMatrix(id, matr);
		}
		
		if(o.altT!=null && totalDepth>IdentityProfile1.writeIsoformDepthThresh){
		o.writeString(id, cc.all_breaks, this.num_sources);
		}
		
	}
	public void process1(CigarCluster cc,   Outputs o,String chrom,int chrom_index){
		String read_count = TranscriptUtils.getString(cc.readCount);
	/*	int startPos, endPos, startPos2, endPos2;
		if(!IdentityProfile1.annotByBreakPosition ){
			startPos = cc.start;
			endPos =cc.end;
			startPos2 =-1;
			endPos2 = -1;
			
		}else{
			startPos = cc.breakSt;
			endPos = cc.breakEnd;
			startPos2 = cc.breakSt2;
			endPos2 = cc.breakEnd2;
		} */
	//	String downstream = annot==null ? null : annot.nextDownstream(endPos, chrom_index);
	//	String upstream = annot==null ? null : annot.nextUpstream(startPos, chrom_index);
	//	String downstream2 = annot==null ||  endPos2<0 ? null : annot.nextDownstream(endPos2, chrom_index);
	//	String upstream2 = annot==null||  startPos2<0  ? null : annot.nextUpstream(startPos2, chrom_index);
	//	if(upstream==null) upstream =  chrom_index+"."+TranscriptUtils.round(startPos, CigarHash2.round)+"";
	//	if(downstream==null) downstream = chrom_index+"."+ TranscriptUtils.round(endPos, CigarHash2.round)+"";
		boolean hasLeaderBreak = TranscriptUtils.coronavirus  ? (cc.breaks.size()>1 &&  annot.isLeader(cc.breaks.get(1)*CigarHash2.round)) : false;

		o.printTranscript(
			cc.id()+"\t"+chrom+"\t"+cc.start+"\t"+cc.end+"\t"+cc.getTypeNme(seqlen)+"\t"+
		//cc.breaks.toString()+"\t"+cc.breaks.hashCode()+"\t"+
	//	cc.breakSt+"\t"+cc.breakEnd+"\t"+cc.breakSt2+"\t"+cc.breakEnd2+"\t"+
		cc.all_breaks.size()+"\t"+cc.numBreaks()+"\t"+(hasLeaderBreak? 1: 0)+"\t"+cc.breaks_hash.secondKey+"\t"+
	//	upstream+"\t"+downstream+"\t"+
	//	upstream2+"\t"+downstream2+"\t"+
		cc.totLen+"\t"+cc.readCountSum()+"\t"+read_count
				+"\t"+cc.getTotDepthSt(true)+"\t"+cc.getTotDepthSt(false)+"\t"+cc.getErrorRatioSt());
		//return TranscriptUtils.round(cc.start, 100)+"."+downstream+"."+upstream+"."+TranscriptUtils.round(seqlen-cc.end, 100);
	}
	
	public void getConsensus(  
			Outputs o, String chrom, int chrom_index
			) throws IOException{
		o.readClusters.close();
		for(Iterator<CigarCluster> it = l.values().iterator(); it.hasNext();) {
			CigarCluster cc = it.next();
			this.process1(cc, o, chrom, chrom_index);
		}
		System.err.println("closing transcripts pw");
	   o.transcriptsP.close();
      
		for(Iterator<CigarCluster> it = l.values().iterator(); it.hasNext();) {
			CigarCluster cc = it.next();
			this.process(cc, o, chrom, chrom_index);
		}
		
		
		
	}

}