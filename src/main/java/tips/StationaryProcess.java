package tips;



public interface StationaryProcess <S> extends Process<S>
{
  public double getStationaryProbability(S state);
}
