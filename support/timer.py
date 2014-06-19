# A with statement based timer.

import timeit

class Timer:
  """
  A with statement based timer.
  Usage:
    t = Timer()
    with t:
      do_something()
    print t.time
  """
    
  def __enter__(self):
    self.start = timeit.default_timer()

  def __exit__(self, type, value, traceback):
    finish = timeit.default_timer()
    self.time = finish - self.start
