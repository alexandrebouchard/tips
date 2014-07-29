package tips;

import java.util.Random;



public interface StationaryProcess <S> extends Process<S>
{
  public double getStationaryProbability(S state);
  public S sampleFromStationary(Random rand);
}
