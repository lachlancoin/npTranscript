package npTranscript.cluster;

import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

import ch.systemsx.cisd.hdf5.IHDF5SimpleWriter;
import japsa.seq.Sequence;
import japsa.seq.SequenceOutputStream;

/**
 * @author Lachlan Coin
 *
 */


public class CigarClusters {
	
	
	
	final double thresh;
	
	CigarClusters(double thresh){
		this.thresh = thresh;
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

	public void getConsensus(Annotation annot, Sequence refseq,  PrintWriter exonP ,
			PrintWriter transcriptsP, SequenceOutputStream seqFasta, IHDF5SimpleWriter  clusterW, 
			
			int[] depth, int num_sources) throws IOException{
		int[] first_last = new int[2];
		int seqlen = refseq.length();
		List<String> str = new ArrayList<String>();
		str.add("pos"); //str.add("ẗype"); 
		for(int j=0; j<num_sources; j++){
			str.add("depth"+j);
		}
		for(int j=0; j<num_sources; j++){
			str.add("errors"+j);
		}
	//	int numcols = str.size();
		//for(int i=0; i<exonP.length; i++){
			exonP.println("ID\tstart\tend\t"+TranscriptUtils.getString("count", num_sources,true));
		//	transcriptsP[i].println("ID,index,start,end,startPos,endPos,totLen,countTotal,"+getString("count", num_sources,true));
			transcriptsP.println("ID\tstart\tend\tbreaks\thash\tstartBreak\tendBreak\tleftGene\trightGene\ttotLen\tcountTotal\t"+TranscriptUtils.getString("count", num_sources,true)
			+"\t"+TranscriptUtils.getString("depth", num_sources, true)+"\t"+TranscriptUtils.getString("errors", num_sources, true)+"\t"+TranscriptUtils.getString("error_ratio", num_sources, true));
			clusterW.writeStringArray("header", str.toArray(new String[0]));
		//}
		
		
		int startPos = 0;
		String sep = "\t";
		int tolerance = 5 ; //for matching to annotation 
		for(Iterator<CigarCluster> it = l.values().iterator(); it.hasNext();) {
			CigarCluster cc = it.next();
			cc.addZeros(); 
		
			 int[][] matr =cc.getClusterDepth(num_sources);
			String id = cc.id+"";
			clusterW.writeIntMatrix(id, matr);
			int[][] exons = cc.getExons( 0.3,10);
		
			String read_count = TranscriptUtils.getString(cc.readCount);
			StringBuffer descline = new StringBuffer();//cc.index+","+read_count);
			StringBuffer subseq= new StringBuffer();
			StringBuffer annotline = new StringBuffer();
			int transcript_len =0;
			
			transcriptsP.println(cc.id+"\t"+cc.start+"\t"+cc.end+"\t"+cc.breaks.toString()+"\t"+cc.breaks.hashCode()+"\t"+
			cc.breakSt+"\t"+cc.breakEnd+"\t"+
			annot.nextUpstream(cc.breakSt, tolerance)+"\t"+annot.nextDownstream(cc.breakEnd, tolerance)+"\t"+
			cc.totLen+"\t"+cc.readCountSum+"\t"+read_count
					+"\t"+cc.getTotDepthSt(true)+"\t"+cc.getTotDepthSt(false)+"\t"+cc.getErrorRatioSt());
			for(int j=0; j<exons.length; j++) {
				int start = exons[j][0];
				int end = exons[j][1];
				exonP.println(id+"\t"+start+"\t"+end+"\t"+read_count);

				
				annotline.append(annot.calcORFOverlap(start, end, first_last, transcript_len));

				int len = end-start+1;
				descline.append(";");
				descline.append(start); descline.append("-"); descline.append(end); descline.append(","); descline.append(len);
				
				//descline.append("|");descline.append(annot.getInfo(first_last[0]));
				//descline.append("|");descline.append(annot.getInfo(first_last[1]));
				subseq.append(refseq.subSequence(start, end).toString());
				
				//System.err.println(subseq.length());
				//System.err.println(subseq);
				//System.err.println("h");
				transcript_len += len;
				//seqline.append(subseq.toString());
			}
			Sequence subseq1 = new Sequence(refseq.alphabet(),subseq.toString().toCharArray(), id);
		//	subseq1.setName(id);
			descline.append(" "); descline.append(annotline);
			subseq1.setDesc(descline.toString());
			subseq1.writeFasta(seqFasta);
		//	seqFasta.println(idline.toString());
		//	seqFasta.println(seqline.toString());
		}
		
	}

}