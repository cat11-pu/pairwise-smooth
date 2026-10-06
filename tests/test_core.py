"""Tests for the curve interpolation and easing kernel (unittest, standard library only)."""

import unittest

from smooth.core import EASING_PRESETS, Bezier, CurveError, ease


class CurveTestCase(unittest.TestCase):
    """Shared comparisons for points and polygons."""

    def assertPointClose(self, first, second, delta=1e-9):
        self.assertAlmostEqual(first[0], second[0], delta=delta)
        self.assertAlmostEqual(first[1], second[1], delta=delta)

    def assertPolygonClose(self, first, second, delta=1e-9):
        self.assertEqual(len(first), len(second))
        for index, point in enumerate(second):
            self.assertPointClose(first[index], point, delta)


class QuadraticTests(CurveTestCase):
    def test_a_quadratic_curve_hits_the_points_between_its_endpoints(self):
        curve = Bezier([(0.0, 0.0), (2.0, 4.0), (4.0, 0.0)])
        self.assertEqual(curve.degree, 2)
        self.assertEqual(curve.bounds(), (0.0, 0.0, 4.0, 4.0))
        self.assertPointClose(curve.point(0.0), (0.0, 0.0))
        self.assertPointClose(curve.point(0.25), (1.0, 1.5))
        self.assertPointClose(curve.point(0.5), (2.0, 2.0))
        self.assertPointClose(curve.point(0.75), (3.0, 1.5))
        self.assertPointClose(curve.point(1.0), (4.0, 0.0))
        previous = -1.0
        for step in range(21):
            column = curve.point(step / 20.0)[0]
            self.assertGreater(column, previous)
            previous = column


class CubicTests(CurveTestCase):
    def test_a_cubic_curve_hits_the_points_between_its_endpoints(self):
        curve = Bezier([(0.0, 0.0), (2.0, 8.0), (6.0, 2.0), (8.0, 0.0)])
        self.assertEqual(curve.degree, 3)
        self.assertPointClose(curve.point(0.0), (0.0, 0.0))
        self.assertPointClose(curve.point(0.25), (1.8125, 3.65625))
        self.assertPointClose(curve.point(0.5), (4.0, 3.75))
        self.assertPointClose(curve.point(0.75), (6.1875, 1.96875))
        self.assertPointClose(curve.point(1.0), (8.0, 0.0))
        self.assertPointClose(curve.velocity(0.0), (6.0, 24.0))
        self.assertPointClose(curve.velocity(1.0), (6.0, -6.0))


class SplitTests(CurveTestCase):
    def test_splitting_a_curve_gives_two_curves_that_meet_on_the_split_point(self):
        cubic = Bezier([(0.0, 0.0), (2.0, 8.0), (6.0, 2.0), (8.0, 0.0)])
        head, tail = cubic.split(0.5)
        self.assertEqual(head.degree, 3)
        self.assertEqual(tail.degree, 3)
        self.assertPolygonClose(
            head.control_points, ((0.0, 0.0), (1.0, 4.0), (2.5, 4.5), (4.0, 3.75))
        )
        self.assertPolygonClose(
            tail.control_points, ((4.0, 3.75), (5.5, 3.0), (7.0, 1.0), (8.0, 0.0))
        )
        self.assertPointClose(head.point(1.0), (4.0, 3.75))
        self.assertPointClose(tail.point(0.0), (4.0, 3.75))

        quadratic = Bezier([(0.0, 0.0), (2.0, 4.0), (4.0, 0.0)])
        head, tail = quadratic.split(0.5)
        self.assertPolygonClose(head.control_points, ((0.0, 0.0), (1.0, 2.0), (2.0, 2.0)))
        self.assertPolygonClose(tail.control_points, ((2.0, 2.0), (3.0, 2.0), (4.0, 0.0)))
        self.assertPointClose(head.point(1.0), (2.0, 2.0))
        self.assertPointClose(tail.point(0.0), (2.0, 2.0))


class TrimTests(CurveTestCase):
    def test_trimming_a_curve_matches_the_same_stretch_of_the_original(self):
        curve = Bezier([(0.0, 0.0), (2.0, 8.0), (6.0, 2.0), (8.0, 0.0)])
        whole = curve.segment(0.0, 1.0)
        self.assertPolygonClose(whole.control_points, curve.control_points)

        middle = curve.segment(0.25, 0.75)
        self.assertPolygonClose(
            middle.control_points,
            ((1.8125, 3.65625), (3.1875, 4.71875), (4.8125, 3.40625), (6.1875, 1.96875)),
        )
        for share in (0.0, 1.0 / 3.0, 0.5, 2.0 / 3.0, 1.0):
            self.assertPointClose(middle.point(share), curve.point(0.25 + 0.5 * share))

        half = curve.segment(0.0, 0.5)
        self.assertPolygonClose(
            half.control_points, ((0.0, 0.0), (1.0, 4.0), (2.5, 4.5), (4.0, 3.75))
        )
        self.assertPointClose(half.point(1.0), curve.point(0.5))

        held = curve.segment(0.75, 0.75)
        for share in (0.0, 0.5, 1.0):
            self.assertPointClose(held.point(share), curve.point(0.75))


class LengthTests(CurveTestCase):
    def test_the_arc_length_grows_with_the_parameter(self):
        line = Bezier([(0.0, 0.0), (3.0, 4.0)])
        self.assertAlmostEqual(line.length(0.0), 0.0, delta=1e-9)
        self.assertAlmostEqual(line.length(0.25), 1.25, delta=1e-9)
        self.assertAlmostEqual(line.length(0.5), 2.5, delta=1e-9)
        self.assertAlmostEqual(line.length(0.75), 3.75, delta=1e-9)
        self.assertAlmostEqual(line.length(1.0), 5.0, delta=1e-9)
        self.assertAlmostEqual(line.length(1.0, steps=8), 5.0, delta=1e-9)
        self.assertAlmostEqual(line.length(1.0, steps=512), 5.0, delta=1e-9)

        curve = Bezier([(0.0, 0.0), (2.0, 8.0), (6.0, 2.0), (8.0, 0.0)])
        self.assertGreater(curve.length(0.25), 0.0)
        self.assertGreater(curve.length(0.5), curve.length(0.25))
        self.assertGreater(curve.length(0.75), curve.length(0.5))
        self.assertGreater(curve.length(1.0), curve.length(0.75))
        self.assertAlmostEqual(
            curve.length(1.0, steps=512), curve.length(1.0, steps=2048), delta=1e-3
        )


class WalkTests(CurveTestCase):
    def test_walking_a_distance_along_the_arc_lands_on_the_matching_point(self):
        line = Bezier([(0.0, 0.0), (3.0, 4.0)])
        parameter, position = line.at_distance(0.0)
        self.assertAlmostEqual(parameter, 0.0, delta=1e-9)
        self.assertPointClose(position, (0.0, 0.0))

        parameter, position = line.at_distance(1.25, steps=8)
        self.assertAlmostEqual(parameter, 0.25, delta=1e-9)
        self.assertPointClose(position, (0.75, 1.0))

        parameter, position = line.at_distance(2.5)
        self.assertAlmostEqual(parameter, 0.5, delta=1e-9)
        self.assertPointClose(position, (1.5, 2.0))

        parameter, position = line.at_distance(5.0, steps=8)
        self.assertAlmostEqual(parameter, 1.0, delta=1e-9)
        self.assertPointClose(position, (3.0, 4.0))

        parameter, position = line.at_distance(6.0, steps=8)
        self.assertAlmostEqual(parameter, 1.0, delta=1e-9)
        self.assertPointClose(position, (3.0, 4.0))


class SampleTests(CurveTestCase):
    def test_evenly_spaced_samples_keep_their_step_along_the_curve(self):
        line = Bezier([(0.0, 0.0), (3.0, 4.0)])
        self.assertPolygonClose(line.samples(2), ((0.0, 0.0), (3.0, 4.0)))
        self.assertPolygonClose(
            line.samples(5),
            ((0.0, 0.0), (0.75, 1.0), (1.5, 2.0), (2.25, 3.0), (3.0, 4.0)),
        )
        self.assertPolygonClose(
            line.samples(9, steps=16),
            tuple((0.375 * index, 0.5 * index) for index in range(9)),
        )

        curve = Bezier([(0.0, 0.0), (2.0, 8.0), (6.0, 8.0), (8.0, 0.0)])
        points = curve.samples(5, steps=512)
        self.assertEqual(len(points), 5)
        self.assertPointClose(points[0], (0.0, 0.0))
        self.assertPointClose(points[-1], (8.0, 0.0))
        self.assertPointClose(points[2], (4.0, 6.0))
        for index in range(5):
            mirrored = (8.0 - points[4 - index][0], points[4 - index][1])
            self.assertPointClose(points[index], mirrored)
        previous = -1.0
        for point in points:
            self.assertGreater(point[0], previous)
            previous = point[0]


class EasingTests(CurveTestCase):
    def test_a_linear_easing_is_the_identity(self):
        for step in range(11):
            self.assertAlmostEqual(ease("linear", step / 10.0), step / 10.0, delta=1e-9)

    def test_easing_presets_keep_the_shape_of_their_timing_curve(self):
        self.assertEqual(
            sorted(EASING_PRESETS), ["ease", "ease-in", "ease-in-out", "ease-out", "linear"]
        )
        self.assertAlmostEqual(ease("ease-in", 0.25), 0.093464651, delta=1e-8)
        self.assertAlmostEqual(ease("ease-in", 0.5), 0.315356813, delta=1e-8)
        self.assertAlmostEqual(ease("ease-in", 0.75), 0.621861869, delta=1e-8)
        self.assertAlmostEqual(ease("ease-out", 0.25), 0.378138131, delta=1e-8)
        self.assertAlmostEqual(ease("ease-out", 0.75), 0.906535349, delta=1e-8)
        self.assertAlmostEqual(ease("ease", 0.5), 0.802403388, delta=1e-8)
        self.assertAlmostEqual(ease("ease-in-out", 0.25), 0.129161931, delta=1e-8)
        self.assertAlmostEqual(ease("ease-in-out", 0.5), 0.5, delta=1e-8)

        for name in sorted(EASING_PRESETS):
            self.assertAlmostEqual(ease(name, 0.0), 0.0, delta=1e-9)
            self.assertAlmostEqual(ease(name, 1.0), 1.0, delta=1e-9)
            previous = -1.0
            for step in range(21):
                value = ease(name, step / 20.0)
                self.assertGreaterEqual(value, previous)
                previous = value


class InputTests(CurveTestCase):
    def test_malformed_input_and_parameters_out_of_range_are_rejected(self):
        with self.assertRaises(CurveError):
            Bezier([(0.0, 0.0)])
        with self.assertRaises(CurveError):
            Bezier([(0.0, 0.0), (1.0, 1.0), (2.0, 2.0), (3.0, 3.0), (4.0, 4.0)])
        with self.assertRaises(CurveError):
            Bezier("(0, 0) (1, 1)")
        with self.assertRaises(CurveError):
            Bezier([(0.0, 0.0, 0.0), (1.0, 1.0)])
        with self.assertRaises(CurveError):
            Bezier([(0.0, 0.0), ("1", 1.0)])
        with self.assertRaises(CurveError):
            Bezier([(0.0, 0.0), (True, 1.0)])

        curve = Bezier([(0.0, 0.0), (3.0, 4.0)])
        for bad in (-0.5, -0.001, 1.001, 2.0):
            with self.assertRaises(CurveError):
                curve.point(bad)
            with self.assertRaises(CurveError):
                curve.length(bad)
        with self.assertRaises(CurveError):
            curve.point("0.5")
        with self.assertRaises(CurveError):
            curve.point(None)
        with self.assertRaises(CurveError):
            curve.split(1.5)
        with self.assertRaises(CurveError):
            curve.segment(0.75, 0.25)
        with self.assertRaises(CurveError):
            curve.segment(-0.25, 0.5)
        with self.assertRaises(CurveError):
            curve.length(1.0, steps=0)
        with self.assertRaises(CurveError):
            curve.length(1.0, steps=1.5)
        with self.assertRaises(CurveError):
            curve.at_distance(-0.5)
        with self.assertRaises(CurveError):
            curve.at_distance(1.0, steps=0)
        with self.assertRaises(CurveError):
            curve.samples(1)
        with self.assertRaises(CurveError):
            curve.samples(4.0)
        with self.assertRaises(CurveError):
            curve.samples(4, steps=0)
        with self.assertRaises(CurveError):
            ease("bounce", 0.5)
        with self.assertRaises(CurveError):
            ease("ease-in", 1.5)
        with self.assertRaises(CurveError):
            ease("ease-in", -0.5)


if __name__ == "__main__":
    unittest.main()
