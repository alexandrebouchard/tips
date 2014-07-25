package tips.bd;

import org.apache.commons.math3.stat.descriptive.SummaryStatistics;
import org.junit.Assert;
import org.junit.Test;

import briefj.collections.Counter;

import tips.TimeIntegratedPathSampler;
import tips.pip.PIPMain;



public class TestCherry
{
  public static void main(String [] args)
  {
    new TestCherry().test();
  }

  @Test
  public void test()
  {
    PIPMain pipMain = new PIPMain();
    pipMain.ensureLinearizationUnique = false;
    pipMain.lambda =  1.3;
    pipMain.mu = 0.5;
    pipMain.bl =  1;
    
    TimeIntegratedPathSampler<Integer> sampler = new TimeIntegratedPathSampler<Integer>(new SimpleBirthDeathPotential(), new ReversibleBDProcess(pipMain.lambda, pipMain.mu));
    
    double num = 0.0, denom = 0.0;
    
    for (int i = 0; i < 1000; i++)
    {
      // sample from PIP    
      pipMain.generateNextData();
      
      // use the len
      final int 
        x = pipMain.getGeneratedEndPoints().sequences().get(PIPMain.ta).length(),
        y = pipMain.getGeneratedEndPoints().sequences().get(PIPMain.tb).length();
      
      // simulate end point problem
      Counter<Integer> rootSamples = sampler.sampleTreeCherry(x, y, 1.0, 1.0);
      rootSamples.normalize();
      denom++;
      for (Integer key : rootSamples.keySet())
        num += rootSamples.getCount(key) * key;
    }
    
    final double analytic =  (pipMain.lambda / pipMain.mu);
    final double mc =  (num/denom);
    Assert.assertEquals(analytic, mc, 0.01);
  }
}
