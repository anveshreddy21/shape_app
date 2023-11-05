import 'dart:math';
import 'dart:ui';
import 'shape_refinement.dart';

Map<String, dynamic> rdpNRWithIndices(List<Offset> points,
    {double alpha = 0.05}) {
  double calculateToleranceBasedOnSize(List<Offset> points, double alpha) {
    double minX = points.reduce((p1, p2) => p1.dx < p2.dx ? p1 : p2).dx;
    double maxX = points.reduce((p1, p2) => p1.dx > p2.dx ? p1 : p2).dx;
    double minY = points.reduce((p1, p2) => p1.dy < p2.dy ? p1 : p2).dy;
    double maxY = points.reduce((p1, p2) => p1.dy > p2.dy ? p1 : p2).dy;

    return alpha *
        sqrt((maxX - minX) * (maxX - minX) + (maxY - minY) * (maxY - minY));
  }

  double pointToLineDistance(Offset p, Offset a, Offset b) {
    if (a == b) {
      return (p - a).distance;
    }

    double A = b.dy - a.dy;
    double B = a.dx - b.dx;
    double C = A * a.dx + B * a.dy;

    return (A * p.dx + B * p.dy - C).abs() / sqrt(A * A + B * B);
  }

  if (points.isEmpty) {
    return {'indices': [], 'points': []};
  }

  Map<int, Offset> simplified = {};
  List<List<int>> stack = [
    [0, points.length - 1]
  ];

  while (stack.isNotEmpty) {
    List<int> segment = stack.removeLast();
    int start = segment[0];
    int end = segment[1];

    double maxDist = 0;
    int index = start;

    for (int i = start + 1; i < end; i++) {
      double dist = pointToLineDistance(points[i], points[start], points[end]);
      if (dist > maxDist) {
        index = i;
        maxDist = dist;
      }
    }

    if (maxDist > calculateToleranceBasedOnSize(points, alpha)) {
      stack.add([index, end]);
      stack.add([start, index]);
    } else {
      simplified[start] = points[start];
      if (stack.isEmpty || stack.last[0] != index) {
        simplified[end] = points[end];
      }
    }
  }

  return {
    'indices': simplified.keys.toList(),
    'points': simplified.values.toList(),
  };
}

//Closure of figure

List<Offset> closure(List<Offset> points) {
  double calculateThresholdBasedOnSize(List<Offset> points,
      {double alpha = 0.1}) {
    double minX = points.reduce((p1, p2) => p1.dx < p2.dx ? p1 : p2).dx;
    double maxX = points.reduce((p1, p2) => p1.dx > p2.dx ? p1 : p2).dx;
    double minY = points.reduce((p1, p2) => p1.dy < p2.dy ? p1 : p2).dy;
    double maxY = points.reduce((p1, p2) => p1.dy > p2.dy ? p1 : p2).dy;

    double diagonal =
        sqrt((maxX - minX) * (maxX - minX) + (maxY - minY) * (maxY - minY));

    return (alpha * diagonal < 10) ? 10 : alpha * diagonal;
  }

  Offset firstPoint = points.first;
  Offset lastPoint = points.last;

  double distance = sqrt(
      (firstPoint.dx - lastPoint.dx) * (firstPoint.dx - lastPoint.dx) +
          (firstPoint.dy - lastPoint.dy) * (firstPoint.dy - lastPoint.dy));

  if (distance < calculateThresholdBasedOnSize(points, alpha: 0.15)) {
    points[points.length - 1] = firstPoint;
    return points;
  }

  return points;
}

//Local angle calculation

List<double> calculateLocalAngles(List<Offset> points, List<int> indices) {
  // Check if the list of indices is too short to process
  if (indices.length <= 1) {
    return []; // Return an empty list
  }

  Offset getVector(Offset p1, Offset p2) {
    return Offset(p2.dx - p1.dx, p2.dy - p1.dy);
  }

  double dotProduct(Offset v1, Offset v2) {
    return v1.dx * v2.dx + v1.dy * v2.dy;
  }

  double magnitude(Offset v) {
    return sqrt(v.dx * v.dx + v.dy * v.dy);
  }

  double angleBetweenVectors(Offset v1, Offset v2) {
    double cosTheta = dotProduct(v1, v2) / (magnitude(v1) * magnitude(v2));
    cosTheta = max(-1, min(1, cosTheta)); // Handle floating point inaccuracies
    return acos(cosTheta) * (180 / pi);
  }

  List<double> angleList = [];

  void computeAngles(int rangeValue) {
    for (int i = 1; i < indices.length - 1; i++) {
      double xAvg1 = 0;
      double xAvg2 = 0;
      double yAvg1 = 0;
      double yAvg2 = 0;
      for (int j = 0; j < rangeValue; j++) {
        xAvg1 += points[indices[i] - j].dx;
        xAvg2 += points[indices[i] + j].dx;
        yAvg1 += points[indices[i] - j].dy;
        yAvg2 += points[indices[i] + j].dy;
      }
      xAvg1 /= rangeValue;
      xAvg2 /= rangeValue;
      yAvg1 /= rangeValue;
      yAvg2 /= rangeValue;

      Offset v1 = getVector(Offset(xAvg1, yAvg1), points[indices[i]]);
      Offset v2 = getVector(points[indices[i]], Offset(xAvg2, yAvg2));
      double angle = angleBetweenVectors(v1, v2);
      angleList.add(angle);
    }
  }

  try {
    computeAngles(5);
  } catch (e) {
    angleList.clear();
    computeAngles(2);
  }

  return angleList;
}

//shape-decider based on the number of RDP points and local angles

Map<String, dynamic> shapeDecider(
    List<Offset> points, List<Offset> rdpPoints, List<double> localAngles) {
  bool isClosed = (rdpPoints.first == rdpPoints.last);
  if (rdpPoints.length == 1) {
    return {'shape': 'Point', 'points': rdpPoints};
  }
  if (rdpPoints.length == 2) {
    return {'shape': 'PolyLine', 'points': rdpPoints};
  }
  if (rdpPoints.length == 3 && !isClosed) {
    return {'shape': 'PolyLine', 'points': rdpPoints};
  }
  if (rdpPoints.length == 4) {
    if (isClosed) {
      List<Offset> triVertices = rdpPoints.sublist(0, 3);
      return {'shape': 'Triangle', 'points': triVertices};
    }
    if (localAngles.every((angle) => angle < 30)) {
      return {'shape': 'Curve', 'points': rdpPoints};
    } else {
      return {'shape': 'PolyLine', 'points': rdpPoints};
    }
  }
  if (rdpPoints.length == 5) {
    if (isClosed && localAngles.every((angle) => angle > 15)) {
      List<Offset> quadVertices = rdpPoints.sublist(0, 4);
      if (isConvexPolygon(quadVertices)) {
        return {'shape': 'Quadrilateral', 'points': quadVertices};
      }
      return {'shape': 'Concave_Polygon', 'points': rdpPoints};
    }
    if (isClosed && localAngles.any((angle) => angle < 15)) {
      return {'shape': 'NoShape', 'points': rdpPoints};
    }
    if (localAngles.every((angle) => angle < 30)) {
      return {'shape': 'Curve', 'points': rdpPoints};
    } else {
      return {'shape': 'PolyLine', 'points': rdpPoints};
    }
  }
  if (rdpPoints.length == 6) {
    if (isClosed && localAngles.every((angle) => angle > 15)) {
      List<Offset> polyVertices = rdpPoints.sublist(0, rdpPoints.length - 1);
      if (isConvexPolygon(polyVertices)) {
        return {'shape': 'Convex_Polygon', 'points': polyVertices};
      }
      return {'shape': 'Concave_Polygon', 'points': polyVertices};
    }
    if (isClosed && localAngles.any((angle) => angle < 15)) {
      return {'shape': 'NoShape', 'points': rdpPoints};
    }
    if (localAngles.every((angle) => angle < 30)) {
      return {'shape': 'Curve', 'points': rdpPoints};
    } else {
      return {'shape': 'PolyLine', 'points': rdpPoints};
    }
  }

  if ([7, 8, 9, 10].contains(rdpPoints.length)) {
    if (isClosed &&
        (localAngles.reduce((a, b) => a + b)) / localAngles.length < 30) {
      Map<String, dynamic> circleResults = leastSquaresCircle(points);
      double h = circleResults['centerX'];
      double k = circleResults['centerY'];
      double r = circleResults['radius'];
      double err = circleResults['error'];

      Map<String, dynamic> minMaxResults = minMaxMethod(points);

      double boxRatio = minMaxResults['boxDim'].dx / minMaxResults['boxDim'].dy;
      if (boxRatio < 1) {
        boxRatio = 1 / boxRatio;
      }

      double aByb = minMaxResults['a'] / minMaxResults['b'];

      if (err < 0.075) {
        return {
          'shape': 'Circle',
          'points': rdpPoints,
          'center': Offset(h, k),
          'centerX': h,
          'centerY': k,
          'radius': r,
          'error': err
        };
      }
      if (err < 0.15 &&
          boxRatio > 0.8 &&
          boxRatio < 1.3 &&
          aByb > 0.8 &&
          aByb < 1.35) {
        return {
          'shape': 'Circle',
          'points': rdpPoints,
          'center': Offset(h, k),
          'centerX': h,
          'centerY': k,
          'radius': r,
          'error': err
        };
      }

      double angle = minMaxResults['tilt'];
      if (angle.abs() <= pi / 12) {
        return {
          'shape': 'Horizontal_Ellipse',
          'points': rdpPoints,
          'center': Offset(h, k),
          'centerX': h,
          'centerY': k,
          'a': minMaxResults['a'],
          'b': minMaxResults['b'],
          'error': err
        };
      } else {
        return {
          'shape': 'Tilted_Ellipse',
          'points': rdpPoints,
          'center': Offset(h, k),
          'centerX': h,
          'centerY': k,
          'a': minMaxResults['a'],
          'b': minMaxResults['b'],
          'angle': angle,
          'error': err
        };
      }
    }
    if (isClosed && localAngles.every((angle) => angle > 15)) {
      List<Offset> polyVertices = rdpPoints.sublist(0, rdpPoints.length - 1);
      if (isConvexPolygon(polyVertices)) {
        return {'shape': 'Convex_Polygon', 'points': polyVertices};
      }
      return {'shape': 'Concave_Polygon', 'points': polyVertices};
    }

    if ((localAngles.reduce((a, b) => a + b)) / localAngles.length < 25) {
      return {'shape': 'Curve', 'points': rdpPoints};
    } else {
      return {'shape': 'PolyLine', 'points': rdpPoints};
    }
  }

  if (rdpPoints.length > 10) {
    if (isClosed && localAngles.every((angle) => angle > 15)) {
      List<Offset> polyVertices = rdpPoints.sublist(0, rdpPoints.length - 1);
      if (isConvexPolygon(polyVertices)) {
        return {'shape': 'Convex_Polygon', 'points': polyVertices};
      }
      return {'shape': 'Concave_Polygon', 'points': polyVertices};
    }
    if ((localAngles.reduce((a, b) => a + b)) / localAngles.length < 25) {
      return {'shape': 'Curve', 'points': rdpPoints};
    } else {
      return {'shape': 'PolyLine', 'points': rdpPoints};
    }
  } else {
    return {'shape': 'NoShape', 'points': rdpPoints};
  }
}

//Circle handling
Map<String, dynamic> leastSquaresCircle(List<Offset> points) {
  int n = points.length;
  if (n == 0) {
    return {'centerX': null, 'centerY': null, 'radius': null, 'error': null};
  }

  // Calculate means
  double meanX =
      (points.map((point) => point.dx).reduce((a, b) => a + b) / n).toDouble();

  double meanY =
      (points.map((point) => point.dy).reduce((a, b) => a + b) / n).toDouble();

  // Calculate terms needed to solve for the parameters
  double suu = points
      .map((point) => pow(point.dx - meanX, 2))
      .reduce((a, b) => a + b)
      .toDouble();

  double suv = points
      .map((point) => (point.dx - meanX) * (point.dy - meanY))
      .reduce((a, b) => a + b)
      .toDouble();
  double svv = points
      .map((point) => pow(point.dy - meanY, 2))
      .reduce((a, b) => a + b)
      .toDouble();
  double suuu = points
      .map((point) => pow(point.dx - meanX, 3))
      .reduce((a, b) => a + b)
      .toDouble();
  double suvv = points
      .map((point) => pow(point.dx - meanX, 2) * (point.dy - meanY))
      .reduce((a, b) => a + b)
      .toDouble();
  double svvv = points
      .map((point) => pow(point.dy - meanY, 3))
      .reduce((a, b) => a + b)
      .toDouble();
  double suvvv = points
      .map((point) => (point.dx - meanX) * pow(point.dy - meanY, 2))
      .reduce((a, b) => a + b)
      .toDouble();
  // Solve the linear system of equations
  double a = 0.5 * (suuu + suvv);
  double b = 0.5 * (svvv + suvvv);
  double uc = (a * svv - b * suv) / (suu * svv - suv * suv);
  double vc = (suu * b - suv * a) / (suu * svv - suv * suv);

  // Convert back from mean-centered coordinates to original coordinates
  double centerX = uc + meanX;
  double centerY = vc + meanY;
  double radius =
      sqrt(pow(centerX - meanX, 2) + pow(centerY - meanY, 2) + (suu + svv) / n);
  double calculateCumulativeDeviation(
      List<Offset> points, double centerX, double centerY, double radius) {
    double totalDeviation = 0;
    for (Offset point in points) {
      double distance =
          sqrt(pow(point.dx - centerX, 2) + pow(point.dy - centerY, 2));
      double deviation = (distance - radius).abs();
      totalDeviation += deviation;
    }
    double cumulativeDeviation = (totalDeviation / points.length).toDouble();
    double error = cumulativeDeviation / radius;
    return error;
  }

  double error = calculateCumulativeDeviation(points, centerX, centerY, radius);
  double minX = points.map((point) => point.dx).reduce((a, b) => min(a, b));
  double maxX = points.map((point) => point.dx).reduce((a, b) => max(a, b));
  double minY = points.map((point) => point.dy).reduce((a, b) => min(a, b));
  double maxY = points.map((point) => point.dy).reduce((a, b) => max(a, b));
  double newCenterX = (minX + maxX) / 2;
  double newCenterY = (minY + maxY) / 2;
  return {
    'centerX': newCenterX,
    'centerY': newCenterY,
    'radius': radius,
    'error': error
  };
}

//Min-Max method for circle and ellipse for better accuracy
Map<String, dynamic> minMaxMethod(List<Offset> points) {
  double distance(Offset p1, Offset p2) {
    double dx = p2.dx - p1.dx;
    double dy = p2.dy - p1.dy;
    return sqrt(dx * dx + dy * dy);
  }

  double cxMin =
      points.reduce((curr, next) => curr.dx < next.dx ? curr : next).dx;
  double cxMax =
      points.reduce((curr, next) => curr.dx > next.dx ? curr : next).dx;
  double cyMin =
      points.reduce((curr, next) => curr.dy < next.dy ? curr : next).dy;
  double cyMax =
      points.reduce((curr, next) => curr.dy > next.dy ? curr : next).dy;
  double cx = (cxMin + cxMax) / 2;
  double cy = (cyMin + cyMax) / 2;
  Offset centre = Offset(cx, cy);
  List<Offset> boundingBox = [
    Offset(cxMin, cyMin),
    Offset(cxMax, cyMin),
    Offset(cxMax, cyMax),
    Offset(cxMin, cyMax)
  ];
  Offset boxDim = Offset(cxMax - cxMin, cyMax - cyMin);
  List<double> dist = points.map((point) => distance(centre, point)).toList();
  double a = dist.reduce(max);
  double b = dist.reduce(min);

  Offset farthestPoint = points[dist.indexOf(a)];
  double tilt;
  if (farthestPoint.dx != centre.dx) {
    tilt =
        atan((farthestPoint.dy - centre.dy) / (farthestPoint.dx - centre.dx));
  } else {
    tilt = pi / 2;
  }

  Offset nearestPoint = points[dist.indexOf(b)];

  return {
    'a': a,
    'b': b,
    'tilt': tilt,
    'boundingBox': boundingBox,
    'boxDim': boxDim,
    'centre': centre,
    'farthestPoint': farthestPoint,
    'nearestPoint': nearestPoint,
  };
}
