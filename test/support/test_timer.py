import os, sys, time, timeit, util
from nose2.tests._common import TestCase

sys.path.append(os.path.join(os.path.dirname(__file__), '..', '..', 'support'))
from timer import Timer, print_time

def do_something():
  time.sleep(0.01)

class TimerTest(TestCase):
  def test_timer(self):
    t = Timer()
    start = timeit.default_timer()
    with t:
      do_something()
    finish = timeit.default_timer()
    time = t.time
    self.assertLessEqual(0.01, time)
    self.assertGreaterEqual(finish - start, time)

  def test_print_time(self):
    stdout = util.CaptureStdout()
    with stdout:
      with print_time('doing', 'something'):
        do_something()
    self.assertEqual(
      'doing something\ndoing something finished in 0.01 second(s)\n',
      stdout.str())
