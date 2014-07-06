package tips.pip;



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

import tips.ImportanceSampler;
import tips.PotPropOptions;
import tips.PotProposal;
import tips.pip.PIPPotential;
import tips.pip.PIPProcess;
import tips.pip.PIPString;




public class TestPIP
{
  
  
  public static void main(String [] args)
  {

    Main pipMain = new Main();
    
    
    for (int repeat = 0; repeat < 10; repeat++)
    {
    
      MSAPoset msa = pipMain.getGeneratedEndPoints();
      
      double exact = Math.exp(exact(pipMain.mu, pipMain.lambda, pipMain.bl, msa));
      
      System.out.println("exact = " + exact);
      
      System.out.println(msa);
      
      ImportanceSampler<PIPString> is = pipMain.buildImportanceSampler();
      is.nParticles = 1;
      SummaryStatistics weightVariance = new SummaryStatistics();
      for (int i = 0; i < 10; i++)
      {
        SummaryStatistics stat = new SummaryStatistics();
        for (int j = 0; j < 100; j++)
        {
          Counter<List<PIPString>> samples = is.sample(pipMain.getStart(), pipMain.getEnd(), pipMain.bl, weightVariance);
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
      n1 = PIPTreeNode.withLabel(Main.ta.toString()),
      n2 = PIPTreeNode.withLabel(Main.tb.toString());
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
    double logPrior = Math.log(pd.probability(msa.sequences().get(Main.ta).length()));
    return logJoint - logPrior;
  }



}
