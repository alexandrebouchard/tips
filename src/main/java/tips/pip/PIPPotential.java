package tips.pip;

import tips.Potential;




public class PIPPotential implements Potential<PIPString>
{

  @Override
  public double get(PIPString proposed, PIPString target)
  {
    int targetNZero = target.zeroes();
    int [] targetPlusses = target.plusses();
    
    int curNPlus = 0;
    int interZeroIdx = 0;
    
    int pot = 0;
    
    for (int curChar : proposed.characters)
    {
      if (curChar == 0)
      {
        pot += Math.abs(targetPlusses[interZeroIdx] - curNPlus);
        
        curNPlus = 0;
        interZeroIdx ++;
      }
      else if (curChar == +1)
        curNPlus++;
      else if (curChar == -1)
        pot ++;
      else 
        throw new RuntimeException();
    }
    pot += Math.abs(targetPlusses[interZeroIdx] - curNPlus);
    
    if (interZeroIdx != targetNZero)
      return Double.POSITIVE_INFINITY;
    
    return pot;
  }
  
}