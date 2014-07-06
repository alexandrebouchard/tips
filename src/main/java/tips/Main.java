package tips;



import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.Set;

import muset.LinearizedAlignment;
import muset.MSAPoset;
import muset.SequenceId;

import org.apache.commons.lang3.tuple.Pair;
import org.apache.commons.math3.distribution.PoissonDistribution;
import org.apache.commons.math3.stat.descriptive.SummaryStatistics;
import org.jgrapht.UndirectedGraph;

import com.google.common.collect.Maps;




import bayonet.graphs.GraphUtils;
import briefj.Indexer;
import briefj.collections.Counter;
import briefj.collections.UnorderedPair;

import tips.pip.PIPPotential;
import tips.pip.PIPProcess;
import tips.pip.PIPString;
import tips.pip.ref.PIPLikelihoodCalculator;
import tips.pip.ref.PIPTreeNode;
import tips.pip.ref.PoissonParameters;




public class Main
{
  private static SequenceId ta = new SequenceId("A"), tb = new SequenceId("B");
  
  public static void main(String [] args)
  {
    double lambda =  2.0;
    double mu = 0.5;
    double bl =  0.3;
    Random genRand = new Random(1);
    
    

    
    PIPProcess process = new PIPProcess(lambda, mu);
    
    for (int repeat = 0; repeat < 10; repeat++)
    {
    
      MSAPoset msa = PIPProcess.keepOnlyEndPts(process.sample(genRand, bl), ta, tb);
      
      double exact = Math.exp(exact(mu, lambda, bl, msa));
      
      System.out.println("exact = " + exact);
      
      System.out.println(msa);
      Pair<PIPString,PIPString> endPoints = getEndPoints(msa, ta, tb);
      PIPPotential pot = new PIPPotential();
      PotProposal<PIPString> prop = new PotProposal<PIPString>(process, pot, new PotPropOptions());
      ImportanceSampler<PIPString> is = new ImportanceSampler<PIPString>(prop, process);
      is.nParticles = 1;
      SummaryStatistics weightVariance = new SummaryStatistics();
      for (int i = 0; i < 10; i++)
      {
        SummaryStatistics stat = new SummaryStatistics();
        for (int j = 0; j < 100; j++)
        {
          Counter<List<PIPString>> samples = is.sample(endPoints.getLeft(), endPoints.getRight(), bl, weightVariance);
          double estimate = is.estimateZ(samples);
          stat.addValue(Math.pow((estimate - exact), 2));
//          System.out.println("approx(" + is.nParticles + ") = " + (estimate));
        }
        System.out.println("mse " + stat.getMean());
        is.nParticles *= 2;
      }
      
      System.out.println();
    }
  }
  
  private static double exact(double mu, double lambda, double bl, MSAPoset msa)
  { 
    LinearizedAlignment linearizedMSA = new LinearizedAlignment(msa);
    PIPTreeNode root = PIPTreeNode.nextUnlabelled();
    PIPTreeNode 
      n1 = PIPTreeNode.withLabel(ta.toString()),
      n2 = PIPTreeNode.withLabel(tb.toString());
    UndirectedGraph<PIPTreeNode, UnorderedPair<PIPTreeNode,PIPTreeNode>> topo = GraphUtils.newUndirectedGraph();
    topo.addVertex(n1);
    topo.addVertex(n2);
    topo.addVertex(root);
    topo.addEdge(n1, root);
    topo.addEdge(n2, root);
    Map<UnorderedPair<PIPTreeNode,PIPTreeNode>, Double> bls = Maps.newLinkedHashMap();
    bls.put(UnorderedPair.of(n1, root), bl/2.0);
    bls.put(UnorderedPair.of(n2, root), bl/2.0);
    
    double [][] trivialRateMtx = new double[][]{{1}};
    Indexer<Character> trivialIndex = new Indexer<Character>();
    trivialIndex.addToIndex(PIPProcess.star);
    
    PoissonParameters poissonParams = new PoissonParameters(trivialIndex, trivialRateMtx, lambda, mu);
    
    PIPLikelihoodCalculator calculator = new PIPLikelihoodCalculator(poissonParams, linearizedMSA, topo, bls, root);
    double logJoint = calculator.computeDataLogProbabilityGivenTree(); 
    
    double statRate = lambda / mu;
    PoissonDistribution pd = new PoissonDistribution(statRate);
    double logPrior = Math.log(pd.probability(msa.sequences().get(ta).length()));
    return logJoint - logPrior;
  }

  public static Pair<PIPString, PIPString> getEndPoints(MSAPoset msa, SequenceId ta,
      SequenceId tb)
  {
    if (msa.sequences().size() != 2)
      throw new RuntimeException();
    Set<Integer> aligned1 = new HashSet<Integer>(), aligned2 = new HashSet<Integer>();
    for (muset.MSAPoset.Column c : msa.columns())
      if (c.getPoints().size() == 2)
      {
        if (c.getPoints().containsKey(ta)) aligned1.add(c.getPoints().get(ta));
        if (c.getPoints().containsKey(tb)) aligned2.add(c.getPoints().get(tb));
      }
    
    List<Integer> cur; PIPString pips1, pips2;
    
    cur = new ArrayList<Integer>(); for (int i = 0; i < msa.sequences().get(ta).length(); i++) cur.add(aligned1.contains(i) ? 0 : -1); pips1 = new PIPString(cur);
    cur = new ArrayList<Integer>(); for (int i = 0; i < msa.sequences().get(tb).length(); i++) cur.add(aligned2.contains(i) ? 0 : +1); pips2 = new PIPString(cur);
    
    return Pair.of(pips1, pips2);
  }

}
