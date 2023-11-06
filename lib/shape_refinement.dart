import 'dart:math';
import 'package:flutter/material.dart';

Offset rotateScalePoint(
    Offset point, Offset pivot, double thetaDegrees, double scale) {
  // Convert angle to radians for Dart
  double theta = thetaDegrees * pi / 180.0;

  // Translate point to the origin:
  double x1 = point.dx - pivot.dx;
  double y1 = point.dy - pivot.dy;

  // Rotate the point around the origin
  double xNew = x1 * cos(theta) - y1 * sin(theta);
  double yNew = x1 * sin(theta) + y1 * cos(theta);

  // Scale the point
  xNew *= scale;
  yNew *= scale;

  // Translate the point back
  xNew += pivot.dx;
  yNew += pivot.dy;

  return Offset(xNew, yNew);
}

double distance(Offset p1, Offset p2) {
  return (p1 - p2).distance;
}

List<Offset> idealTriangle(List<Offset> vertices, double thresholdAngle) {
  Offset p1 = vertices[0], p2 = vertices[1], p3 = vertices[2];

  double calculateSlope(double x1, double y1, double x2, double y2) {
    if (x1 == x2) return double.infinity;
    return (y2 - y1) / (x2 - x1);
  }

  double findRotationAngle(List<Offset> vertices) {
    var A = vertices[0], B = vertices[1], C = vertices[2];
    List slopes = [
      [calculateSlope(A.dx, A.dy, B.dx, B.dy), 'AB'],
      [calculateSlope(B.dx, B.dy, C.dx, C.dy), 'BC'],
      [calculateSlope(C.dx, C.dy, A.dx, A.dy), 'CA'],
    ];

    double smallestAngle = double.infinity;
    double rotationAngle = 0.0;
    for (var item in slopes) {
      double slope = item[0];
      double angleToHorizontal;

      if (slope == double.infinity) {
        angleToHorizontal = 90.0;
      } else {
        angleToHorizontal = atan(slope).abs() * (180 / pi);
      }

      if (angleToHorizontal < smallestAngle) {
        smallestAngle = angleToHorizontal;
        rotationAngle = atan(slope) * (180 / pi);
      }
    }
    return rotationAngle;
  }

  double rotationAngle = findRotationAngle(vertices);
  if (rotationAngle.abs() > thresholdAngle) {
    return vertices;
  }

  double dx0 = p2.dx - p1.dx;
  double dy0 = p2.dy - p1.dy;
  double dx1 = p3.dx - p1.dx;
  double dy1 = p3.dy - p1.dy;
  double crossProduct = dx0 * dy1 - dy0 * dx1;
  if (crossProduct < 0) {
    Offset temp = p2;
    p2 = p3;
    p3 = temp;
  }

  Offset p2New = rotateScalePoint(p2, p1, -rotationAngle, 1.0);
  Offset p3New = rotateScalePoint(p3, p1, -rotationAngle, 1.0);
  return [p1, p2New, p3New];
}

List<double> triangleAngles(Offset A, Offset B, Offset C) {
  double a = distance(B, C);
  double b = distance(A, C);
  double c = distance(A, B);

  double alpha = acos((b * b + c * c - a * a) / (2 * b * c));
  double beta = acos((a * a + c * c - b * b) / (2 * a * c));
  double gamma = acos((a * a + b * b - c * c) / (2 * a * b));

  return [alpha, beta, gamma].map((rad) => rad * (180 / pi)).toList();
}

String classifyTriangle(List<double> angles) {
  if (angles.every((angle) => angle >= 45 && angle <= 75)) {
    return "Equilateral Triangle";
  }

  for (var angle in angles) {
    if (angle >= 80 && angle <= 100) {
      var otherAngles = angles.where((a) => a != angle).toList();
      if ((otherAngles[0] - otherAngles[1]).abs() <= 5) {
        return "Isosceles Right-Angled Triangle";
      }
      return "Right-Angled Triangle";
    }
  }

  if ((angles[0] - angles[1]).abs() <= 8 ||
      (angles[1] - angles[2]).abs() <= 8 ||
      (angles[2] - angles[0]).abs() <= 8) {
    return "Isosceles Triangle";
  }

  return "Scalene Triangle";
}

List<Offset> idealScaleneTriangle(List<Offset> vertices, List<double> angles) {
  return idealTriangle(vertices, 8.0);
}

List<Offset> idealEquilateralTriangle(
    List<Offset> vertices, List<double> angles) {
  Offset midpoint(Offset p1, Offset p2) {
    return Offset((p1.dx + p2.dx) / 2, (p1.dy + p2.dy) / 2);
  }

  double factor = (distance(vertices[0], vertices[1]) +
          distance(vertices[1], vertices[2]) +
          distance(vertices[2], vertices[0])) /
      3;
  List<Offset> equilateralVertices(
      Offset baseStart, Offset baseEnd, Offset thirdPoint) {
    Offset baseEnd2 = rotateScalePoint(
        baseEnd, baseStart, 0, factor / distance(baseStart, baseEnd));
    double dx = baseEnd2.dx - baseStart.dx;
    double dy = baseEnd2.dy - baseStart.dy;
    double L = sqrt(pow(dx, 2) + pow(dy, 2));

    // Compute height of equilateral triangle
    double height = L * sqrt(3) / 2;

    // Compute direction to position the third vertex
    int direction = ((baseEnd.dx - baseStart.dx) *
                    (thirdPoint.dy - midpoint(baseStart, baseEnd).dy) -
                (baseEnd.dy - baseStart.dy) *
                    (thirdPoint.dx - midpoint(baseStart, baseEnd).dx)) >
            0
        ? 1
        : -1;

    Offset thirdVertex = Offset(
      midpoint(baseStart, baseEnd2).dx - direction * dy * height / L,
      midpoint(baseStart, baseEnd2).dy + direction * dx * height / L,
    );

    List<Offset> verticesNew = [baseStart, baseEnd2, thirdVertex];
    return verticesNew;
  }

  double horizontalScore(Offset point1, Offset point2) {
    // Compute a score indicating how horizontal a line segment is.
    return (point2.dx - point1.dx).abs();
  }

  // Find the most horizontal side
  List<double> scores = List.generate(
      3, (i) => horizontalScore(vertices[i], vertices[(i + 1) % 3]));
  int baseIndex = scores.indexOf(scores.reduce(max));
  Offset baseStart = vertices[baseIndex];
  Offset baseEnd = vertices[(baseIndex + 1) % 3];
  Offset thirdPoint = vertices[(baseIndex + 2) % 3];

  // Generate ideal equilateral triangle with the found base
  List<Offset> idealVertices1 =
      equilateralVertices(baseStart, baseEnd, thirdPoint);

  return idealTriangle(idealVertices1, 8.0);
}

List<Offset> idealRightIsoscelesTriangle(
    List<Offset> vertices, List<double> angles) {
  double maxAngle = angles.reduce(max);
  int maxAngleIndex = angles.indexOf(maxAngle);
  Offset p = vertices[maxAngleIndex];
  Offset a = vertices[(maxAngleIndex + 1) % 3];
  Offset b = vertices[(maxAngleIndex + 2) % 3];
  double scale = distance(p, b) / distance(p, a);
  double theta = maxAngle - 90;
  Offset aNew = theta != 0 ? rotateScalePoint(a, p, theta, scale) : a;
  return idealTriangle(
      [p, aNew, b], 8.0); // idealTriangle function should be implemented
}

List<Offset> idealRightTriangle(List<Offset> vertices, List<double> angles) {
  double maxAngle = angles.reduce(max);
  int maxAngleIndex = angles.indexOf(maxAngle);
  Offset p = vertices[maxAngleIndex];
  Offset a = vertices[(maxAngleIndex + 1) % 3];
  Offset b = vertices[(maxAngleIndex + 2) % 3];
  double crossProduct = (a.dx - p.dx) * (b.dy - p.dy) -
      (a.dy - p.dy) * (b.dx - p.dx); // cross product of vectors ap and bp
  if (crossProduct < 0) {
    Offset temp = a;
    a = b;
    b = temp;
  }
  double theta = maxAngle - 90;
  Offset aNew = rotateScalePoint(a, p, theta, 1.0);
  return idealTriangle(
      [p, aNew, b], 8.0); // idealTriangle function should be implemented
}

//Isosceles case
List<Offset> idealIsoscelesTriangle(
    List<Offset> vertices, List<double> angles) {
  List<double> angleDif = List.generate(
    3,
    (i) => (angles[(i + 1) % 3] - angles[(i + 2) % 3]).abs(),
  );
  int pivotIndex = angleDif.indexOf(angleDif.reduce(min));
  Offset pivot = vertices[pivotIndex];
  Offset a = vertices[(pivotIndex + 1) % 3];
  Offset b = vertices[(pivotIndex + 2) % 3];
  double pa = distance(pivot, a);
  double pb = distance(pivot, b);
  double length = (pa + pb) / 2;
  double scaleA = length / pa;
  double scaleB = length / pb;

  // Ensure A is to the left of B when viewed from pivot
  double crossProduct = (a.dx - pivot.dx) * (b.dy - pivot.dy) -
      (a.dy - pivot.dy) * (b.dx - pivot.dx);
  if (crossProduct < 0) {
    Offset temp = a;
    a = b;
    b = temp;
    double tempScale = scaleA;
    scaleA = scaleB;
    scaleB = tempScale;
  }

  Offset aNew = rotateScalePoint(a, pivot, 0, scaleA);
  double alphaPlusBeta =
      (angles[(pivotIndex + 1) % 3] + angles[(pivotIndex + 2) % 3]);
  Offset bNew = rotateScalePoint(aNew, pivot, 180 - alphaPlusBeta, scaleB);
  return idealTriangle(
      [pivot, aNew, bNew], 8.0); // idealTriangle function should be implemented
}

//function that takes 3 points and returns the vertices of the final triangle depending on the triangle type
List<Offset> finalTriangle(List<Offset> vertices) {
  List<double> angleTriangle =
      triangleAngles(vertices[0], vertices[1], vertices[2]);
  String triangleType = classifyTriangle(angleTriangle);
  List<Offset> idealTriangle;
  switch (triangleType) {
    case "Equilateral Triangle":
      idealTriangle = idealEquilateralTriangle(vertices, angleTriangle);

      break;
    case "Isosceles Right-Angled Triangle":
      var idealTriangleIso = idealIsoscelesTriangle(vertices, angleTriangle);
      idealTriangle = idealRightTriangle(
          idealTriangleIso,
          triangleAngles(
              idealTriangleIso[0], idealTriangleIso[1], idealTriangleIso[2]));

      break;
    case "Right-Angled Triangle":
      idealTriangle = idealRightTriangle(vertices, angleTriangle);

      // print(triangleAngles(idealTriangle[0], idealTriangle[1], idealTriangle[2])); // Debug print
      break;
    case "Isosceles Triangle":
      idealTriangle = idealIsoscelesTriangle(vertices, angleTriangle);

      // print(triangleAngles(idealTriangle[0], idealTriangle[1], idealTriangle[2])); // Debug print
      break;
    case "Scalene Triangle":
      idealTriangle = idealScaleneTriangle(vertices, angleTriangle);

      break;
    default:
      idealTriangle = vertices;
      break;
  }
  return idealTriangle;
}

//Refinement of quadrilaterals

//Quad angles
List<double> quadAngles(List<Offset> vertices) {
  if (vertices.length != 4) {
    throw ArgumentError("Input should contain 4 vertices.");
  }

  List<double> vector(Offset p1, Offset p2) {
    return [p2.dx - p1.dx, p2.dy - p1.dy];
  }

  double dotProduct(List<double> v1, List<double> v2) {
    return v1[0] * v2[0] + v1[1] * v2[1];
  }

  double magnitude(List<double> v) {
    return sqrt(pow(v[0], 2) + pow(v[1], 2));
  }

  double angleBetweenVectors(List<double> v1, List<double> v2) {
    double cosineAngle = dotProduct(v1, v2) / (magnitude(v1) * magnitude(v2));
    // Ensure the cosine does not fall out of the range due to floating point errors.
    cosineAngle = cosineAngle.clamp(-1.0, 1.0);
    double angle = acos(cosineAngle);
    return angle * (180 / pi);
  }

  List<double> angles = [];
  for (int i = 0; i < 4; i++) {
    List<double> v1 = vector(vertices[i], vertices[(i + 1) % 4]);
    List<double> v2 = vector(vertices[i], vertices[(i + 3) % 4]);
    angles.add(angleBetweenVectors(v1, v2));
  }

  return angles;
}

//Quad type
String quadType(List<Offset> vertices, List<double> angles) {
  bool quadIsConvex(List<Offset> vertices) {
    bool ccw(Offset A, Offset B, Offset C) {
      return (C.dy - A.dy) * (B.dx - A.dx) > (B.dy - A.dy) * (C.dx - A.dx);
    }

    bool intersect(Offset A, Offset B, Offset C, Offset D) {
      return ccw(A, C, D) != ccw(B, C, D) && ccw(A, B, C) != ccw(A, B, D);
    }

    return !intersect(vertices[0], vertices[1], vertices[2], vertices[3]);
  }

  if (!quadIsConvex(vertices)) {
    return "concave_quadrilateral";
  }

  double distance(Offset p1, Offset p2) {
    return sqrt(pow(p2.dx - p1.dx, 2) + pow(p2.dy - p1.dy, 2));
  }

  List<double> sides =
      List.generate(4, (i) => distance(vertices[i], vertices[(i + 1) % 4]));
  double meanSide = sides.reduce((a, b) => a + b) / 4;
  bool sideCriteria =
      sides.every((side) => (side - meanSide).abs() / meanSide <= 0.12);
  bool angleCriteria = angles.every((angle) => (angle - 90).abs() <= 7);

  if (sideCriteria && angleCriteria) {
    return "square";
  }
  if (sideCriteria) {
    return "rhombus";
  }
  if (angleCriteria) {
    return "rectangle";
  }
  if ((angles[0] - angles[2]).abs() <= 10 &&
      (angles[1] - angles[3]).abs() <= 10) {
    return "parallelogram";
  }
  for (int i = 0; i < 4; i++) {
    if (175 <= angles[i] + angles[(i + 1) % 4] &&
        angles[i] + angles[(i + 1) % 4] <= 185 &&
        175 <= angles[(i + 2) % 4] + angles[(i + 3) % 4] &&
        angles[(i + 2) % 4] + angles[(i + 3) % 4] <= 185) {
      return "trapezium";
    }
  }
  return "irregular_quadrilateral";
}

//Ideal Quadrilateral
List<Offset> idealQuadrilateral(List<Offset> vertices, double thresholdAngle) {
  Offset p1 = vertices[0], p2 = vertices[1], p3 = vertices[2], p4 = vertices[3];

  double calculateSlope(Offset p1, Offset p2) {
    if (p1.dx == p2.dx) {
      return double.infinity;
    }
    return (p2.dy - p1.dy) / (p2.dx - p1.dx);
  }

  double findRotationAngle(List<Offset> vertices) {
    var slopes = [
      calculateSlope(vertices[0], vertices[1]),
      calculateSlope(vertices[1], vertices[2]),
      calculateSlope(vertices[2], vertices[3]),
      calculateSlope(vertices[3], vertices[0]),
    ];

    double smallestAngle = double.infinity;
    double rotationAngle = 0;
    for (var slope in slopes) {
      double angleToHorizontal;
      if (slope == double.infinity) {
        angleToHorizontal = 90;
      } else {
        angleToHorizontal = atan(slope).abs() * 180 / pi;
      }

      if (angleToHorizontal < smallestAngle) {
        smallestAngle = angleToHorizontal;
        rotationAngle = atan(slope) * 180 / pi;
      }
    }

    return rotationAngle;
  }

  double rotAngle = findRotationAngle(vertices);
  if (rotAngle.abs() > thresholdAngle) {
    return vertices;
  }

  Offset p2New = rotateScalePoint(p2, p1, -rotAngle, 1);
  Offset p3New = rotateScalePoint(p3, p1, -rotAngle, 1);
  Offset p4New = rotateScalePoint(p4, p1, -rotAngle, 1);
  return [p1, p2New, p3New, p4New];
}

//Ideal Square
List<Offset> perfectSquare(List<Offset> vertices) {
  int ccw(Offset A, Offset B, Offset C) {
    return (C.dy - A.dy) * (B.dx - A.dx) > (B.dy - A.dy) * (C.dx - A.dx)
        ? 1
        : -1;
  }

  Offset vertex4 = rotateScalePoint(vertices[1], vertices[0],
      ccw(vertices[0], vertices[1], vertices[2]) * 90.0, 1);
  Offset vertex3 = Offset((vertices[1].dx + vertex4.dx - vertices[0].dx),
      (vertices[1].dy + vertex4.dy - vertices[0].dy));

  return [vertices[0], vertices[1], vertex3, vertex4];
}

//Ideal Rhombus
List<Offset> perfectRhombus(List<Offset> vertices) {
  Offset midpoint(Offset p1, Offset p2) {
    return Offset((p1.dx + p2.dx) / 2, (p1.dy + p2.dy) / 2);
  }

  int ccw(Offset A, Offset B, Offset C) {
    return (C.dy - A.dy) * (B.dx - A.dx) > (B.dy - A.dy) * (C.dx - A.dx)
        ? 1
        : -1;
  }

  Offset A = vertices[0], B = vertices[1], C = vertices[2];

  // Calculate midpoint of AC
  Offset M = midpoint(A, C);

  Offset BPrime = rotateScalePoint(
      A, M, ccw(M, A, B) * 90.0, distance(M, B) / distance(M, A));

  // Reflect B' over M to get D'
  Offset DPrime = Offset(M.dx - (BPrime.dx - M.dx), M.dy - (BPrime.dy - M.dy));

  return [A, BPrime, C, DPrime];
}

//Ideal Rectangle
List<Offset> perfectRectangle(List<Offset> vertices) {
  Offset? intersection(Offset p1, Offset p2, Offset q1, Offset q2) {
    double A1 = p2.dy - p1.dy;
    double B1 = p1.dx - p2.dx;
    double C1 = A1 * p1.dx + B1 * p1.dy;

    double A2 = q2.dy - q1.dy;
    double B2 = q1.dx - q2.dx;
    double C2 = A2 * q1.dx + B2 * q1.dy;

    double det = A1 * B2 - A2 * B1;
    if (det == 0) {
      // Lines are parallel, no intersection
      return null;
    }

    double x = (B2 * C1 - B1 * C2) / det;
    double y = (A1 * C2 - A2 * C1) / det;
    return Offset(x, y);
  }

  Offset A = vertices[0], B = vertices[1], C = vertices[2], D = vertices[3];

  // Calculate intersection of AC and BD
  Offset? O = intersection(A, C, B, D);
  if (O == null) {
    throw Exception(
        'The diagonals AC and BD are parallel and do not intersect.');
  }

  // Find the longest distance from O
  double OAPrime = vertices.map((vertex) => distance(O, vertex)).reduce(max);

  Offset APrime = rotateScalePoint(A, O, 0, OAPrime / distance(O, A));

  // Calculate C' using A' and center O
  Offset CPrime = Offset(2 * O.dx - APrime.dx, 2 * O.dy - APrime.dy);

  // Scale B about O to get B'
  Offset BPrime = rotateScalePoint(B, O, 0, OAPrime / distance(O, B));

  // Reflect B' over O to get D'
  Offset DPrime = Offset(2 * O.dx - BPrime.dx, 2 * O.dy - BPrime.dy);

  return [APrime, BPrime, CPrime, DPrime];
}

//Ideal Parallelogram
List<Offset> perfectParallelogram(List<Offset> vertices) {
  Offset midpoint(Offset p1, Offset p2) {
    return Offset((p1.dx + p2.dx) / 2, (p1.dy + p2.dy) / 2);
  }

  Offset A = vertices[0];
  Offset B = vertices[1];
  Offset C = vertices[2];

  Offset M = midpoint(A, C);
  Offset DPrime = Offset(2 * M.dx - B.dx, 2 * M.dy - B.dy);

  return [A, B, C, DPrime];
}

//Ideal Trapezium
List<Offset> perfectTrapezium(List<Offset> vertices) {
  double calculateSlope(Offset p1, Offset p2) {
    if (p1.dx == p2.dx) {
      return double.infinity;
    }
    return (p2.dy - p1.dy) / (p2.dx - p1.dx);
  }

  List<double> slopes = [];
  List<double> angles = [];

  for (int i = 0; i < vertices.length; i++) {
    slopes
        .add(calculateSlope(vertices[i], vertices[(i + 1) % vertices.length]));
    angles.add(atan(slopes.last));
  }

  Offset A = vertices[0];
  Offset B = vertices[1];
  Offset C = vertices[2];
  Offset D = vertices[3];

  if ((angles[0].abs() - angles[2].abs()).abs() >
      (angles[1].abs() - angles[3].abs()).abs()) {
    A = vertices[3];
    B = vertices[0];
    C = vertices[1];
    D = vertices[2];
  }

  Offset DPrime = Offset(
    C.dx - (B.dx - A.dx) * (distance(C, D) / distance(B, A)),
    C.dy - (B.dy - A.dy) * (distance(C, D) / distance(B, A)),
  );

  return [A, B, C, DPrime];
}

//Final quadrilateral
List<Offset> finalQuadrilateral(List<Offset> vertices) {
  List<double> angles = quadAngles(vertices);
  String quad = quadType(vertices, angles);
  List<Offset> idealQuad;
  switch (quad) {
    case "square":
      idealQuad = perfectSquare(vertices);
      break;
    case "rhombus":
      idealQuad = perfectRhombus(vertices);
      break;
    case "rectangle":
      idealQuad = perfectRectangle(vertices);
      break;
    case "parallelogram":
      idealQuad = perfectParallelogram(vertices);
      break;
    case "trapezium":
      idealQuad = perfectTrapezium(vertices);
      break;
    default:
      idealQuad = vertices;
      break;
  }
  return idealQuadrilateral(idealQuad, 8.0);
}

//Polygon convex or not
bool isConvexPolygon(List<Offset> vertices) {
  if (vertices.length < 4) return true; // Triangles are always convex.

  bool sign = false;
  for (int i = 0; i < vertices.length; i++) {
    Offset A = vertices[i];
    Offset B = vertices[(i + 1) % vertices.length];
    Offset C = vertices[(i + 2) % vertices.length];

    Offset AB = B - A;
    Offset BC = C - B;

    double crossProductZ = AB.dx * BC.dy - AB.dy * BC.dx;
    if (i == 0) {
      // For the first time, save the sign of z of the cross product
      sign = crossProductZ > 0;
    } else {
      // If the sign changes, then polygon is concave.
      if ((crossProductZ > 0) != sign) return false;
    }
  }

  return true; // No sign change found, polygon is convex.
}

//Irregular to Regular polygon handling
List<Offset> finalPolygon(List<Offset> vertices) {
  final int n = vertices.length; // Number of vertices
  if (n < 3) return vertices; // Not enough vertices to form a polygon
  
  // Calculate the centroid of the polygon
  final centroid = vertices.fold<Offset>(
        Offset.zero,
        (prev, element) => prev + element,
      ) /
      n.toDouble();

  // Calculate all side lengths
  List<double> sideLengths = [];
  for (int i = 0; i < n; i++) {
    sideLengths.add((vertices[i] - vertices[(i + 1) % n]).distance);
  }

  // Determine if all sides are approximately equal
  final double maxSide = sideLengths.reduce(max);
  final double minSide = sideLengths.reduce(min);
  const double toleranceRatio = 1.8; // Tolerance of 80%
  if (maxSide / minSide > toleranceRatio) {
    // Not a regular polygon, return the original vertices
    return vertices;
  }
  // Calculate the angle of each vertex for a regular polygon
  final double angleIncrement = 2 * pi / n;
  List<Offset> regularVertices = [];

  // Align first vertex with positive x-axis
  double startAngle = atan2(
    vertices.first.dy - centroid.dy,
    vertices.first.dx - centroid.dx,
  );
  double d = 0.5 * maxSide + 0.5 * minSide;
  double r = d / (2 * sin(pi / n));

  // Generate the vertices of a regular polygon
  for (int i = 0; i < n; i++) {
    double angle = startAngle + angleIncrement * i;
    regularVertices.add(Offset(
      centroid.dx + r * cos(angle),
      centroid.dy + r * sin(angle),
    ));
  }

  return alignPolygon(regularVertices);
}

List<Offset> alignPolygon(List<Offset> vertices) {
  Offset centroid(List<Offset> vertices) {
    double centerX = 0.0, centerY = 0.0;
    for (var vertex in vertices) {
      centerX += vertex.dx;
      centerY += vertex.dy;
    }
    return Offset(centerX / vertices.length, centerY / vertices.length);
  }

  double angleWithHorizontal(Offset a, Offset b) {
    return atan2(b.dy - a.dy, b.dx - a.dx);
  }

  double determineRotationAngle(List<double> angles) {
    double minDistance = double.infinity;
    double rotationAngle = 0.0;

    for (var angle in angles) {
      // Normalize the angle between -π to π
      double normalizedAngle = angle.remainder(2 * pi);
      if (normalizedAngle > pi) {
        normalizedAngle -= 2 * pi;
      } else if (normalizedAngle < -pi) {
        normalizedAngle += 2 * pi;
      }

      // Check distance to horizontal axis (0 or π)
      double distanceToHorizontal =
          min(normalizedAngle.abs(), (normalizedAngle - pi).abs());

      // Check distance to vertical axis (-π/2 or π/2)
      double distanceToVertical = min(
          (normalizedAngle - pi / 2).abs(), (normalizedAngle + pi / 2).abs());

      // Determine if the current angle is closer to the horizontal or vertical axis
      double currentMinDistance, currentRotationAngle;
      if (distanceToHorizontal <= distanceToVertical) {
        currentMinDistance = distanceToHorizontal;
        currentRotationAngle = (normalizedAngle.abs() <= pi / 2)
            ? -normalizedAngle
            : (normalizedAngle > 0)
                ? pi - normalizedAngle
                : -pi - normalizedAngle;
      } else {
        currentMinDistance = distanceToVertical;
        currentRotationAngle = (normalizedAngle > 0)
            ? pi / 2 - normalizedAngle
            : -pi / 2 - normalizedAngle;
      }

      // Check if the current angle provides a closer alignment than what we have so far
      if (currentMinDistance < minDistance) {
        minDistance = currentMinDistance;
        rotationAngle = currentRotationAngle;
      }
    }

    return rotationAngle;
  }

  List<Offset> rotatePolygon(
      List<Offset> vertices, double angle, Offset center) {
    return vertices.map((vertex) {
      final dx = vertex.dx - center.dx;
      final dy = vertex.dy - center.dy;
      final newX = center.dx + dx * cos(angle) - dy * sin(angle);
      final newY = center.dy + dx * sin(angle) + dy * cos(angle);
      return Offset(newX, newY);
    }).toList();
  }

  final int vertexCount = vertices.length;
  final Offset center = centroid(vertices);
  List<double> principalDiagonalAngles = [];

  if (vertexCount.isEven) {
    for (int i = 0; i < vertexCount / 2; i++) {
      principalDiagonalAngles.add(angleWithHorizontal(
          vertices[i], vertices[(i + vertexCount ~/ 2) % vertexCount]));
    }
  } else {
    for (int i = 0; i < vertexCount; i++) {
      final oppositeMidIndex = (i + vertexCount ~/ 2) % vertexCount;
      final oppositeMidPoint = Offset(
        (vertices[(oppositeMidIndex + 1) % vertexCount].dx +
                vertices[oppositeMidIndex].dx) /
            2,
        (vertices[(oppositeMidIndex + 1) % vertexCount].dy +
                vertices[oppositeMidIndex].dy) /
            2,
      );
      principalDiagonalAngles
          .add(angleWithHorizontal(vertices[i], oppositeMidPoint));
    }
  }
  
  final rotationAngle = determineRotationAngle(principalDiagonalAngles);
  return rotatePolygon(vertices, rotationAngle, center);
}

/*
//Not using this right now
//Perfecting the regular polygon

List<Offset> alignBaseWithHorizontal(List<Offset> vertices,
    {double threshold = 0.28}) {
  // Sort the vertices based on the 'y' value to get the bottom ones.
  List<Offset> bottomVertices = List<Offset>.from(vertices)
    ..sort((a, b) => b.dy.compareTo(a.dy));
  // Take the bottom two points which form the potential base.
  
  Offset p1 = bottomVertices[0];
  Offset p2 = bottomVertices[1];
  
  // Calculate the angle of the base with the horizontal axis.
  double angleWithHorizontal = atan2(p2.dy - p1.dy, p2.dx - p1.dx);

  double rotationAngle = -angleWithHorizontal;
  

  if (angleWithHorizontal.abs() > threshold) {
    // The base is already aligned with the horizontal axis.
    return vertices;
  }

  // Calculate the angle needed to rotate the base to align with the horizontal axis.

  // Find the centroid of the polygon to use as the pivot for rotation.
  Offset centroid = vertices.fold(
        Offset.zero,
        (Offset sum, Offset v) => sum + v,
      ) /
      vertices.length.toDouble();

  // Rotate all vertices around the centroid by the calculated angle.
  List<Offset> rotatedVertices = vertices.map((vertex) {
    double dx = vertex.dx - centroid.dx;
    double dy = vertex.dy - centroid.dy;

    return Offset(
      centroid.dx + dx * cos(rotationAngle) - dy * sin(rotationAngle),
      centroid.dy + dx * sin(rotationAngle) + dy * cos(rotationAngle),
    );
  }).toList();

  return rotatedVertices;
}

*/

//Curve Handling

Offset cmr(Offset p0, Offset p1, Offset p2, Offset p3, double t) {
  double t2 = t * t;
  double t3 = t2 * t;
  List<double> a1 = [];
  List<double> a2 = [];
  List<double> a3 = [];
  List<double> a4 = [];

  a1.add(-0.5 * p0.dx + 1.5 * p1.dx - 1.5 * p2.dx + 0.5 * p3.dx);
  a2.add(p0.dx - 2.5 * p1.dx + 2 * p2.dx - 0.5 * p3.dx);
  a3.add(-0.5 * p0.dx + 0.5 * p2.dx);
  a4.add(p1.dx);

  a1.add(-0.5 * p0.dy + 1.5 * p1.dy - 1.5 * p2.dy + 0.5 * p3.dy);
  a2.add(p0.dy - 2.5 * p1.dy + 2 * p2.dy - 0.5 * p3.dy);
  a3.add(-0.5 * p0.dy + 0.5 * p2.dy);
  a4.add(p1.dy);

  double x = a1[0] * t3 + a2[0] * t2 + a3[0] * t + a4[0];
  double y = a1[1] * t3 + a2[1] * t2 + a3[1] * t + a4[1];

  return Offset(x, y);
}
