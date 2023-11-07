import 'package:flutter/material.dart';
import 'analysis.dart';
import 'shape_refinement.dart';

void main() => runApp(MyApp());

class MyApp extends StatelessWidget {
  @override
  Widget build(BuildContext context) {
    return MaterialApp(
      title: 'Drawing App',
      theme: ThemeData(
        primarySwatch: Colors.blue,
      ),
      home: MyHomePage(),
    );
  }
}

class MyHomePage extends StatefulWidget {
  @override
  _MyHomePageState createState() => _MyHomePageState();
}

class _MyHomePageState extends State<MyHomePage> {
  List<Offset> currentStroke = [];
  List<Offset> rdpPoints = [];
  List<int> rdpIndices = [];
  List<double> localAngles = [];
  Map<String, dynamic> shape = {};
  bool showOriginalDrawing = true; // To toggle visibility
  double strokeWidth = 5.0; // Default stroke width

  Offset _getLocalPosition(Offset globalPosition) {
    final RenderBox renderBox = context.findRenderObject() as RenderBox;
    Offset localPosition = renderBox.globalToLocal(globalPosition);
    // Adjust the localPosition by the AppBar height
    localPosition = localPosition.translate(0, -AppBar().preferredSize.height);
    return localPosition;
  }

  @override
  Widget build(BuildContext context) {
    return Scaffold(
      appBar: AppBar(
        title: Text('Drawing App'),
      ),
      body: SafeArea(
        child: Column(
          children: <Widget>[
            Expanded(
              child: GestureDetector(
                onPanUpdate: (details) {
                  setState(() {
                    Offset localPosition =
                        _getLocalPosition(details.globalPosition);
                    currentStroke.add(localPosition);
                  });
                },
                onPanDown: (details) {
                  setState(() {
                    currentStroke.clear();
                    rdpPoints.clear();
                    rdpIndices.clear();
                    localAngles.clear();
                    shape.clear();
                    Offset localPosition =
                        _getLocalPosition(details.globalPosition);
                    currentStroke.add(localPosition);
                  });
                },
                onPanEnd: (details) {
                  setState(() {
                    var result = rdpNRWithIndices(currentStroke, alpha: 0.05);
                    rdpPoints = result['points'] as List<Offset>;
                    rdpIndices = result['indices'] as List<int>;
                    localAngles =
                        calculateLocalAngles(currentStroke, rdpIndices);
                    rdpPoints = closure(currentStroke, rdpPoints);
                    shape = shapeDecider(currentStroke, rdpPoints, localAngles);
                  });
                },
                child: CustomPaint(
                  painter: MyPainter(
                    showOriginalDrawing ? currentStroke : [],
                    showOriginalDrawing ? rdpPoints : [],
                    shape,
                    strokeWidth,
                    showOriginalDrawing,
                  ),
                  size: Size.infinite,
                ),
              ),
            ),
            // Slider and Checkbox
            Padding(
              padding: const EdgeInsets.all(8.0),
              child: Row(
                mainAxisAlignment: MainAxisAlignment.spaceEvenly,
                children: <Widget>[
                  Expanded(
                    child: Slider(
                      value: strokeWidth,
                      min: 1.0,
                      max: 10.0,
                      onChanged: (double value) {
                        setState(() {
                          strokeWidth = value;
                        });
                      },
                    ),
                  ),
                  Checkbox(
                    value: showOriginalDrawing,
                    onChanged: (bool? value) {
                      setState(() {
                        showOriginalDrawing = value ?? true;
                      });
                    },
                  ),
                  Text('Show Raw Drawing'),
                ],
              ),
            ),
            // Shape Type Display
            Padding(
              padding: const EdgeInsets.symmetric(vertical: 8.0),
              child: Text(
                'Shape: ${shape['shape'] ?? 'NoShape'}',
                style: Theme.of(context).textTheme.titleMedium,
              ),
            ),
          ],
        ),
      ),
    );
  }
}

class MyPainter extends CustomPainter {
  final List<Offset> points;
  final List<Offset> rdpPoints;
  final Map<String, dynamic> shape;
  final double strokeWidth;
  final bool showOriginalDrawing;

  MyPainter(this.points, this.rdpPoints, this.shape, this.strokeWidth,
      this.showOriginalDrawing);

  @override
  void paint(Canvas canvas, Size size) {
    Paint paint = Paint()
      ..color = Colors.black
      ..strokeCap = StrokeCap.round
      ..strokeWidth = strokeWidth;

    if (showOriginalDrawing || shape['shape'] == 'NoShape') {
      for (int i = 0; i < points.length - 1; i++) {
        canvas.drawLine(points[i], points[i + 1], paint);
      }
    }

    Paint rdpPaint = Paint()
      ..color = Colors.red
      ..strokeCap = StrokeCap.round
      ..strokeWidth = 7.0;

    for (var point in rdpPoints) {
      canvas.drawCircle(point, 6.0, rdpPaint);
    }

    Paint shapePaint = Paint()
      ..style = PaintingStyle.stroke
      ..color = Colors.blue
      ..strokeCap = StrokeCap.round
      ..strokeWidth = 6.0;

    if (shape['shape'] == 'Circle') {
      canvas.drawCircle(shape['center'], shape['radius'], shapePaint);
    } else if (shape['shape'] == 'Horizontal_Ellipse') {
      double a = shape['a'];
      a = 2 * a;
      double b = shape['b'];
      b = 2 * b;
      canvas.drawOval(
          Rect.fromCenter(center: shape['center'], width: a, height: b),
          shapePaint);
    } else if (shape['shape'] == 'Tilted_Ellipse') {
      double a = shape['a']; // Semi-major axis
      double b = shape['b']; // Semi-minor axis
      Offset center = shape['center']; // Center of the ellipse
      double angle = shape['angle']; // Rotation angle in radians

      // Create a rect where the ellipse will be drawn
      Rect rect = Rect.fromCenter(center: center, width: a * 2, height: b * 2);

      // Save the current canvas state
      canvas.save();

      // Move the canvas to the center of where you want your ellipse
      canvas.translate(center.dx, center.dy);

      // Apply the rotation
      canvas.rotate(angle);

      // Move the canvas back to the original position
      canvas.translate(-center.dx, -center.dy);

      // Draw the oval
      canvas.drawOval(rect, shapePaint);

      // Restore the canvas to its original state
      canvas.restore();
    } else if (shape['shape'] == 'Triangle') {
      List<Offset> vertices = finalTriangle(shape['points']);
      for (int i = 0; i < 3; i++) {
        canvas.drawLine(vertices[i], vertices[(i + 1) % 3], shapePaint);
      }
    } else if (shape['shape'] == 'Quadrilateral') {
      List<Offset> vertices = finalQuadrilateral(shape['points']);
      for (int i = 0; i < 4; i++) {
        canvas.drawLine(vertices[i], vertices[(i + 1) % 4], shapePaint);
      }
    } else if (shape['shape'] == 'PolyLine') {
      List<Offset> vertices = shape['points'];

      for (int i = 0; i < vertices.length - 1; i++) {
        canvas.drawLine(
            vertices[i], vertices[(i + 1) % vertices.length], shapePaint);
      }
    } else if (shape['shape'] == 'Convex_Polygon') {
      List<Offset> polyVertices = finalPolygon(shape['points']);
      for (int i = 0; i < polyVertices.length; i++) {
        canvas.drawLine(polyVertices[i],
            polyVertices[(i + 1) % polyVertices.length], shapePaint);
      }
    } else if (shape['shape'] == 'Concave_Polygon') {
      List<Offset> vertices = shape['points'];

      for (int i = 0; i < vertices.length; i++) {
        canvas.drawLine(
            vertices[i], vertices[(i + 1) % vertices.length], shapePaint);
      }
    } else if (shape['shape'] == 'Curve') {
      List<Offset> controlPoints = shape['points'];

      final spline = CatmullRomSpline(controlPoints, tension: 0);

      final path = Path();
      path.moveTo(controlPoints.first.dx, controlPoints.first.dy);
      for (double t = 0.0; t < 1.0; t += 0.01) {
        final point = spline.transform(t);
        path.lineTo(point.dx, point.dy);
      }
      path.moveTo(controlPoints.last.dx, controlPoints.last.dy);

      // Draw the path on the canvas.
      canvas.drawPath(path, shapePaint);
    }
  }

  @override
  bool shouldRepaint(covariant CustomPainter oldDelegate) => true;
}
