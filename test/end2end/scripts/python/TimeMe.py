import gc
import timeit
import time


class TimeMe:
    """Small class to time portions of the code.
     It has two main usages:
     1. In a with block, of the type:
       t = TimeMe()
       with t:
          doMyThings()
       print("Took {} seconds".format(t.interval)
     2. With multiple ticks and tocs:
       t = TimeMe()
       t.tic() # stores the marker 'start'
       ...
       t.tic("halfway") # stores the marker 'halfway'
       ...
       secondsElapsed = t.toc() # get the seconds from start
       secondsElapsedFromHalfway= t.toc() # get the seconds from halfway
    """

    def __init__(self, timer=None, disable_gc=False, verbose=False):
        if timer is None:
            timer = timeit.default_timer
        self.timer = timer
        self.disable_gc = disable_gc
        self.verbose = verbose
        self.start = self.end = self.interval = None
        self.markers = dict()

    def __enter__(self):
        if self.disable_gc:
            self.gc_state = gc.isenabled()
            gc.disable()
        self.start = self.timer()
        return self

    def __exit__(self, *args):
        self.end = self.timer()
        if self.disable_gc and self.gc_state:
            gc.enable()
        self.interval = self.end - self.start
        if self.verbose:
            print('time taken: %f seconds' % self.interval)

    def print(self, end='\n'):
        """ Print the time elapsed in the with block """
        print("%.4fs" % self.interval, end=end, flush=True)

    def tocStr(self, markerName: str = "start", restart: bool = False):
        """ Get a string representation of the time elapsed from the specified marker """
        return "%.4fs" % self.toc(markerName)

    def tick(self, markerName: str = "start"):
        self.markers[markerName] = self.timer()

    def toc(self, markerName: str = "start", restart: bool = False):
        return self.timer() - self.markers[markerName]
