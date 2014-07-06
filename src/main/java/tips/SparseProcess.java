package tips;

import briefj.collections.Counter;




public interface SparseProcess<S>
{
  public Counter<S> rates(S point);
}