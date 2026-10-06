"""Bezier interpolation and easing curves over plain floats.

A curve is a control polygon of two, three or four points.  It is evaluated
with the Bernstein weights of its degree, split and trimmed through the de
Casteljau construction, and measured with a polyline arc length that can also
be inverted, so a curve can be walked and resampled at even steps.  The easing
presets are unit cubic Bezier timing curves whose x coordinate is inverted by
bisection, the way CSS timing functions are read.

Everything is arithmetic on floats: no clock, no randomness, no I/O.
"""

from math import hypot

EASING_PRESETS = {
    "linear": (0.0, 0.0, 1.0, 1.0),
    "ease": (0.25, 0.1, 0.25, 1.0),
    "ease-in": (0.42, 0.0, 1.0, 1.0),
    "ease-out": (0.0, 0.0, 0.58, 1.0),
    "ease-in-out": (0.42, 0.0, 0.58, 1.0),
}

DEFAULT_STEPS = 64
SOLVER_STEPS = 60


class CurveError(ValueError):
    """Raised for malformed control points, bad parameters or unknown names."""


def _check_number(value, what):
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise CurveError("%s must be a number, got %r" % (what, value))
    return float(value)


def _check_pair(value, what):
    if not isinstance(value, (list, tuple)) or len(value) != 2:
        raise CurveError("%s must be a pair of numbers, got %r" % (what, value))
    return _check_number(value[0], what), _check_number(value[1], what)


def _check_unit(value, what):
    parameter = _check_number(value, what)
    if parameter < 0.0 or parameter > 1.0:
        raise CurveError("%s must lie between 0 and 1, got %r" % (what, value))
    return parameter


def _check_steps(value):
    if isinstance(value, bool) or not isinstance(value, int):
        raise CurveError("a step count must be an integer, got %r" % (value,))
    if value < 1:
        raise CurveError("a step count must be positive, got %r" % (value,))
    return value


def _check_count(value):
    if isinstance(value, bool) or not isinstance(value, int):
        raise CurveError("a sample count must be an integer, got %r" % (value,))
    if value < 2:
        raise CurveError("a sample count must be at least two, got %r" % (value,))
    return value


def _lerp(first, second, weight):
    return (
        first[0] + (second[0] - first[0]) * weight,
        first[1] + (second[1] - first[1]) * weight,
    )


def _spread(first, second):
    return hypot(second[0] - first[0], second[1] - first[1])


def _de_casteljau(points, weight):
    """The successive levels of the de Casteljau construction at one parameter."""
    levels = [tuple(points)]
    current = tuple(points)
    while len(current) > 1:
        current = tuple(
            _lerp(current[index], current[index + 1], weight)
            for index in range(len(current) - 1)
        )
        levels.append(current)
    return levels


class Bezier(object):
    """A Bezier curve of degree one, two or three over the unit parameter range."""

    __slots__ = ("_points",)

    def __init__(self, control_points):
        if not isinstance(control_points, (list, tuple)):
            raise CurveError("control points must be a sequence, got %r" % (control_points,))
        points = tuple(_check_pair(point, "a control point") for point in control_points)
        if len(points) < 2:
            raise CurveError("a curve needs at least two control points, got %d" % (len(points),))
        if len(points) > 4:
            raise CurveError("at most four control points are supported, got %d" % (len(points),))
        self._points = points

    def __repr__(self):
        return "Bezier(%r)" % (self._points,)

    @property
    def control_points(self):
        """The control polygon of the curve."""
        return self._points

    @property
    def degree(self):
        """The degree of the curve."""
        return len(self._points) - 1

    def bounds(self):
        """The bounding box of the control polygon, as (left, top, right, bottom)."""
        columns = [point[0] for point in self._points]
        rows = [point[1] for point in self._points]
        return min(columns), min(rows), max(columns), max(rows)

    def point(self, parameter):
        """The position on the curve at the given parameter."""
        t = _check_unit(parameter, "a curve parameter")
        first, last = self._points[0], self._points[-1]
        if self.degree == 1:
            return _lerp(first, last, t)
        if self.degree == 2:
            middle = self._points[1]
            near = (1.0 - t) * (1.0 - t)
            hump = (1.0 - t) * t
            far = t * t
            return (
                near * first[0] + hump * middle[0] + far * last[0],
                near * first[1] + hump * middle[1] + far * last[1],
            )
        second, third = self._points[1], self._points[2]
        near = (1.0 - t) ** 3
        rise = 3.0 * (1.0 - t) * t * t
        hump = 3.0 * (1.0 - t) * (1.0 - t) * t
        far = t ** 3
        return (
            near * first[0] + rise * second[0] + hump * third[0] + far * last[0],
            near * first[1] + rise * second[1] + hump * third[1] + far * last[1],
        )

    def velocity(self, parameter):
        """The first derivative of the curve at the given parameter."""
        t = _check_unit(parameter, "a curve parameter")
        points = self._points
        if self.degree == 1:
            return points[1][0] - points[0][0], points[1][1] - points[0][1]
        if self.degree == 2:
            first, middle, last = points
            return (
                2.0 * ((1.0 - t) * (middle[0] - first[0]) + t * (last[0] - middle[0])),
                2.0 * ((1.0 - t) * (middle[1] - first[1]) + t * (last[1] - middle[1])),
            )
        first, second, third, last = points
        left = second[0] - first[0], second[1] - first[1]
        middle = third[0] - second[0], third[1] - second[1]
        right = last[0] - third[0], last[1] - third[1]
        near = (1.0 - t) * (1.0 - t)
        far = t * t
        return (
            3.0 * (near * left[0] + 2.0 * (1.0 - t) * t * middle[0] + far * right[0]),
            3.0 * (near * left[1] + 2.0 * (1.0 - t) * t * middle[1] + far * right[1]),
        )

    def split(self, parameter):
        """The two curves that meet at the given parameter and cover the same arc."""
        t = _check_unit(parameter, "a split parameter")
        levels = _de_casteljau(self._points, t)
        left = Bezier([level[0] for level in levels])
        right = Bezier([level[-1] for level in levels])
        return left, right

    def segment(self, start, end):
        """The piece of the curve between two parameters, as a curve of its own."""
        first = _check_unit(start, "a trim parameter")
        last = _check_unit(end, "a trim parameter")
        if last < first:
            raise CurveError("a trim interval must not run backwards, got %r" % ((first, last),))
        if first == last:
            held = _de_casteljau(self._points, first)[-1][0]
            return Bezier([held] * len(self._points))
        levels = _de_casteljau(self._points, first)
        tail = [level[-1] for level in reversed(levels)]
        span = last - first
        levels = _de_casteljau(tail, span)
        return Bezier([level[0] for level in levels])

    def length(self, limit=1.0, steps=DEFAULT_STEPS):
        """The polyline length of the curve up to the given parameter."""
        limit = _check_unit(limit, "an arc length limit")
        steps = _check_steps(steps)
        total = 0.0
        previous = self._points[0]
        for index in range(1, steps + 1):
            current = self.point(index / steps)
            total += _spread(previous, current)
            previous = current
        return total

    def _table(self, steps):
        """Cumulative polyline lengths, one entry per parameter step."""
        totals = [0.0]
        previous = self._points[0]
        for index in range(1, steps + 1):
            current = self.point(index / steps)
            totals.append(totals[-1] + _spread(previous, current))
            previous = current
        return totals

    def at_distance(self, distance, steps=DEFAULT_STEPS):
        """The parameter and position reached after walking a given arc length."""
        walked = _check_number(distance, "an arc length")
        if walked < 0.0:
            raise CurveError("an arc length must not be negative, got %r" % (distance,))
        steps = _check_steps(steps)
        totals = self._table(steps)
        if walked >= totals[-1]:
            parameter = (steps - 1.0) / steps
            return parameter, self.point(parameter)
        index = 0
        while totals[index + 1] < walked:
            index += 1
        span = totals[index + 1] - totals[index]
        share = 0.0 if span <= 0.0 else (walked - totals[index]) / span
        parameter = (index + share) / steps
        return parameter, self.point(parameter)

    def samples(self, count, steps=DEFAULT_STEPS):
        """Points spread evenly along the arc length of the curve."""
        wanted = _check_count(count)
        steps = _check_steps(steps)
        total = self.length(1.0, steps)
        points = []
        for index in range(wanted):
            distance = total * index / wanted
            points.append(self.at_distance(distance, steps)[1])
        return tuple(points)


def _timing_axis(first, second, u):
    """One axis of a unit cubic timing curve with the given two inner ordinates."""
    return (
        3.0 * (1.0 - u) * (1.0 - u) * u * first
        + 3.0 * (1.0 - u) * u * u * second
        + u * u * u
    )


def _timing_parameter(x1, y1, x2, y2, target):
    """The curve parameter whose x coordinate reaches the target."""
    low, high = 0.0, 1.0
    for _ in range(SOLVER_STEPS):
        middle = (low + high) / 2.0
        if _timing_axis(y1, y2, middle) < target:
            low = middle
        else:
            high = middle
    return (low + high) / 2.0


def _timing_points(name):
    if name not in EASING_PRESETS:
        raise CurveError(
            "unknown easing name %r, expected one of %r"
            % (name, tuple(sorted(EASING_PRESETS)))
        )
    return EASING_PRESETS[name]


def ease(name, progress):
    """The eased progress of a named timing curve at progress in [0, 1]."""
    t = _check_unit(progress, "an eased progress")
    x1, y1, x2, y2 = _timing_points(name)
    u = _timing_parameter(x1, y1, x2, y2, t)
    return _timing_axis(y1, y2, u)
