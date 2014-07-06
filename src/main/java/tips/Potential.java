package tips;



public interface Potential<S>
{
  public double get(S proposed, S target);
}
