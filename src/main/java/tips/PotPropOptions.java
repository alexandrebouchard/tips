package tips;

import java.util.Random;

import briefj.opt.Option;



public class PotPropOptions
{
  @Option public boolean automatic = true;
  @Option public double stopPr = 0.95;
  @Option public double greed = 2.0/3.0;
  @Option public String specialSymbol = PotProposal.SPECIAL_SYMBOL;

  public double randPr(Random rand)
  {
    double rInt = rand.nextInt(4) + 2;
    return 1.0 - Math.pow(0.5, rInt);
  }
}