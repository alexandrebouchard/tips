package tips.pip;

import java.util.List;

import muset.MSAPoset;
import muset.MSAPoset.Column;

import org.apache.commons.math3.stat.descriptive.SummaryStatistics;
import org.junit.Assert;
import org.junit.Test;

import com.google.common.collect.Lists;



public class TestExpectedNChanges
{
  public static void main(String [] args)
  {
    new TestExpectedNChanges().test();
  }

  @Test 
  public void test()
  {
    PIPMain pipMain = new PIPMain();
    pipMain.ensureLinearizationUnique = false;
    pipMain.lambda =  1.3;
    pipMain.mu = 0.5;
    pipMain.bl =  1;
    
    SummaryStatistics 
      nChangesStatistics = new SummaryStatistics(),
      statioLength = new SummaryStatistics(),
      nDelsStats = new SummaryStatistics(),
      nInsStats = new SummaryStatistics(),
      nGhostsStats = new SummaryStatistics();
    for (int i = 0; i < 100000; i++)
    {
      pipMain.generateNextData();
      
      MSAPoset fullGeneratedPath = pipMain.getFullGeneratedPath();
//      System.out.println(fullGeneratedPath);
      int nChanges = nChanges(fullGeneratedPath);
      nChangesStatistics.addValue(nChanges);
      statioLength.addValue(pipMain.getGeneratedEndPoints().sequences().get(PIPMain.tb).length());
      
      int nDels = nDels(pipMain.getGeneratedEndPoints());
      nDelsStats.addValue(nDels);
      
      int nIns = nIns(pipMain.getGeneratedEndPoints());
      nInsStats.addValue(nIns);
      
      int nGhosts = 
        pipMain.getFullGeneratedPath().columns().size() - 
        pipMain.getGeneratedEndPoints().columns().size();
      nGhostsStats.addValue(nGhosts);
      
      if (2*nGhosts + nIns + nDels != nChanges)
        throw new RuntimeException();
    }
    
    Assert.assertEquals(formula(pipMain.lambda, pipMain.mu), nChangesStatistics.getMean(), 0.01);
    
//    System.out.println("nChange stats");
//    System.out.println(nChangesStatistics);
//    System.out.println("Analytic: " + formula(pipMain.lambda, pipMain.mu));
//    
//    System.out.println("---");
//    System.out.println("statioLen stats");
//    System.out.println(statioLength);
//    System.out.println("Analytic: " + (pipMain.lambda/pipMain.mu));
//    
//    System.out.println("---");
//    System.out.println("nDels");
//    System.out.println(nDelsStats);
//    System.out.println("Analytic: " + ((pipMain.lambda/pipMain.mu) * (1.0 - Math.exp(-pipMain.mu))));
//    
//    System.out.println("---");
//    System.out.println("nIns");
//    System.out.println(nInsStats);
//    System.out.println("Analytic: " + ((pipMain.lambda/pipMain.mu) * (1.0 - Math.exp(-pipMain.mu))));
//    
//    System.out.println("---");
//    System.out.println("nGhosts");
//    System.out.println(nGhostsStats);
//    System.out.println("Analytic: " +  (pipMain.lambda * (1.0 - (1.0 - Math.exp(pipMain.mu))/pipMain.mu)));
//    
//    System.out.println("----");
//    System.out.println("sanityCheck");
//    System.out.println(nChangesStatistics.getMean());
//    System.out.println(nInsStats.getMean() + nDelsStats.getMean() + 2.0 * nGhostsStats.getMean());
  }

  private int nIns(MSAPoset generatedEndPoints)
  {
    int result = 0;
    for (Column c : generatedEndPoints.columns())
    {
      if (c.getPoints().containsKey(PIPMain.tb) && !c.getPoints().containsKey(PIPMain.ta))
        result++;
    }
    return result;
  }

  private int nDels(MSAPoset generatedEndPoints)
  {
    int result = 0;
    for (Column c : generatedEndPoints.columns())
    {
      if (c.getPoints().containsKey(PIPMain.ta) && !c.getPoints().containsKey(PIPMain.tb))
        result++;
    }
    return result;
  }

  private double formula(double lambda, double mu)
  {
    return 2.0 * lambda;
  }

  private int nChanges(MSAPoset fullGeneratedPath)
  {
    if (fullGeneratedPath.nSequences() == 2)
    {
      List<String> seqs = Lists.newArrayList(fullGeneratedPath.sequences().values());
      return seqs.get(0).length() == seqs.get(1).length() ? 0 : 1;
    }
    return (fullGeneratedPath.nSequences() - 1);
  }
}
