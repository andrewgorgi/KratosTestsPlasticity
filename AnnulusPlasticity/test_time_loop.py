"""Regression: dt=0.00005 must produce 40 steps, not 41."""
import unittest
from types import SimpleNamespace
from annulus_study import Analysis


class TimeLoopTests(unittest.TestCase):
    def test_accumulated_time_stops_at_final_step(self):
        for dt,expected_steps in ((.00005,40),(.0001,20),(.0002,10)):
            state = SimpleNamespace(time=0.,end_time=.002)
            steps = 0
            while Analysis.KeepAdvancingSolutionLoop(state):
                state.time += dt
                steps += 1
                self.assertLessEqual(steps,expected_steps)
            self.assertEqual(steps,expected_steps)
            self.assertAlmostEqual(state.time,state.end_time,places=15)


if __name__ == '__main__':
    unittest.main()
