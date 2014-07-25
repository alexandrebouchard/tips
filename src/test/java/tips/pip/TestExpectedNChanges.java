package tips.pip;

import java.util.List;

import muset.MSAPoset;
import muset.MSAPoset.Column;

import org.apache.commons.math3.stat.descriptive.SummaryStatistics;
import org.junit.Assert;
import org.junit.Test;

import com.google.common.collect.Lists;


/**
 * Check that the expected number of changes (inserts + deletes) in 
 * a branch of length 1 is 2*lambda
 * 
 * @author Alexandre Bouchard (alexandre.bouchard@gmail.com)
 */
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
      statioLength = new SummaryStatistics();
    for (int i = 0; i < 100000; i++)
    {
      pipMain.generateNextData();
      
      MSAPoset fullGeneratedPath = pipMain.getFullGeneratedPath();
      int nChanges = nChanges(fullGeneratedPath);
      nChangesStatistics.addValue(nChanges);
      statioLength.addValue(pipMain.getGeneratedEndPoints().sequences().get(PIPMain.tb).length());
      
      int nDels = nDels(pipMain.getGeneratedEndPoints());
      
      int nIns = nIns(pipMain.getGeneratedEndPoints());
      
      int nGhosts = 
        pipMain.getFullGeneratedPath().columns().size() - 
        pipMain.getGeneratedEndPoints().columns().size();
      
      if (2*nGhosts + nIns + nDels != nChanges)
        throw new RuntimeException();
    }
    
    Assert.assertEquals(formula(pipMain.lambda, pipMain.mu), nChangesStatistics.getMean(), 0.01);
    Assert.assertEquals(pipMain.lambda / pipMain.mu, statioLength.getMean(), 0.01);
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
