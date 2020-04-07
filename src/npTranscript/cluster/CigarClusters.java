package npTranscript.cluster;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

import ch.systemsx.cisd.hdf5.HDF5Factory;
import ch.systemsx.cisd.hdf5.IHDF5SimpleWriter;
import japsa.seq.Sequence;
import japsa.seq.SequenceOutputStream;

/**
 * @author Lachlan Coin
 *
 */


public class CigarClusters {
	
	
	
//	final double thresh;
	
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


	
	public int matchCluster(CigarCluster c1,  int source_index, int num_sources) {
	
		int clusterID=-1;
	
		
		if(!l.containsKey(c1.breaks)){
			CigarCluster newc = new CigarCluster(l.keySet().size(), num_sources, c1, source_index);
			
			clusterID = newc.id;
			l.put(newc.breaks, newc);
			//System.err.println("new cluster " +" "+newc.id+" "+index);

		}else{
			CigarCluster clust = l.get(c1.breaks);
			clust.merge(c1);
			clusterID = clust.id;
		}
		return clusterID;
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
	public String process(CigarCluster cc,   PrintWriter transcriptsP, SequenceOutputStream seqFasta1,IHDF5SimpleWriter  clusterW ,
			Map<String, List<CigarHash>> geneToHash, String altID){
		cc.addZeros(seqlen); 
		
		 int[][] matr =cc.getClusterDepth(num_sources);
		String id = cc.id+"";
		clusterW.writeIntMatrix(id, matr);
		
		String read_count = TranscriptUtils.getString(cc.readCount);
		String downstream = annot.nextDownstream(cc.breakEnd);
		String upstream = annot.nextUpstream(cc.breakSt);
		String ccid = altID==null ? cc.id+"": altID;
		transcriptsP.println(
			ccid+"\t"+cc.start+"\t"+cc.end+"\t"+cc.getTypeNme(seqlen)+"\t"+
		cc.breaks.toString()+"\t"+cc.breaks.hashCode()+"\t"+
		cc.breakSt+"\t"+cc.breakEnd+"\t"+
		upstream+"\t"+downstream+"\t"+
		cc.totLen+"\t"+cc.readCountSum+"\t"+read_count
				+"\t"+cc.getTotDepthSt(true)+"\t"+cc.getTotDepthSt(false)+"\t"+cc.getErrorRatioSt());
		return TranscriptUtils.round(cc.start, 100)+"."+downstream+"."+upstream+"."+TranscriptUtils.round(seqlen-cc.end, 100);
	}
	
	public void getConsensus(  
			PrintWriter transcriptsP,  SequenceOutputStream seqFasta1, File  outfile2
			) throws IOException{
		List<String> str = new ArrayList<String>();
		str.add("pos"); //str.add("ẗype"); 
		for(int j=0; j<num_sources; j++){
			str.add("depth"+j);
		}
		for(int j=0; j<num_sources; j++){
			str.add("errors"+j);
		}
		//	exonP.println("ID\tstart\tend\t"+TranscriptUtils.getString("count", num_sources,true));
			String transcriptP_header = "ID\tstart\tend\ttype_nme\tbreaks\thash\tstartBreak\tendBreak\tleftGene\trightGene\ttotLen\tcountTotal\t"+TranscriptUtils.getString("count", num_sources,true)
			+"\t"+TranscriptUtils.getString("depth", num_sources, true)+"\t"+TranscriptUtils.getString("errors", num_sources, true)
			+"\t"+TranscriptUtils.getString("error_ratio", num_sources, true);
			
			transcriptsP.println(transcriptP_header);
		//	transcriptsP1.println(transcriptP_header);
			IHDF5SimpleWriter clusterW = 
					 HDF5Factory.open(outfile2);
			clusterW.writeStringArray("header", str.toArray(new String[0]));
        Map<String, List<CigarHash>> geneToHash = new HashMap<String, List<CigarHash>>();
		for(Iterator<CigarCluster> it = l.values().iterator(); it.hasNext();) {
			CigarCluster cc = it.next();
			String key = this.process(cc, transcriptsP, seqFasta1, clusterW, geneToHash, null);
			if(geneToHash!=null){
				List<CigarHash> l;
				if(!geneToHash.containsKey(key)) geneToHash.put(key, l = new ArrayList<CigarHash>());
				else  l = geneToHash.get(key);
				l.add(cc.breaks);
			}
		}
		clusterW.close();
		/*{
			IHDF5SimpleWriter clusterW1 = 
					 HDF5Factory.open(outfile2_1);
			clusterW1.writeStringArray("header", str.toArray(new String[0]));
			List<CigarCluster> l1 = new ArrayList<CigarCluster>();
			for(Iterator<List<CigarHash>> it = geneToHash.values().iterator(); it.hasNext();){
				List<CigarHash> cl = it.next();
				CigarCluster cc = l.get(cl.get(0));
				StringBuffer sb = new StringBuffer();
				sb.append(cc.id);
				for(int i=1; i<cl.size(); i++){
					CigarCluster cc2 = l.remove(cl.get(i));
					sb.append(";"+cc2.id);
					cc.merge(cc2);
				}
			
				this.process(cc, transcriptsP1, null, clusterW1, null, sb.toString());
				//Integer.
			}
			clusterW1.close();
			
		}*/
		
	}

}