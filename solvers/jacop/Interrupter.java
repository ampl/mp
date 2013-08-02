import JaCoP.search.ConsistencyListener;

class InterruptSearch extends RuntimeException {}

/**
 * ConsistencyListener that allows to interrupt search.
 */
public class Interrupter implements ConsistencyListener {
  private ConsistencyListener[] consistencyListeners;
  private long data;

  private static native boolean stop(long data);
  
  public Interrupter(long data) {
    this.data = data;
  }

  public void setChildrenListeners(ConsistencyListener[] children) {
    consistencyListeners = children;
  }

  public void setChildrenListeners(ConsistencyListener child) {
    consistencyListeners = new ConsistencyListener[] {child};
  }

  public boolean executeAfterConsistency(boolean consistent) {
    if (stop(data))
      throw new InterruptSearch();
    if (consistencyListeners != null) {
      boolean code = false;
      for (int i = 0; i < consistencyListeners.length; i++)
        code |= consistencyListeners[i].executeAfterConsistency(consistent);
      consistent = code;
    }
    return consistent;
  }
}
