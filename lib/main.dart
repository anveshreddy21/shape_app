import 'dart:math';

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
                'Shape: ${shape['shape'] ?? ''}',
                style: Theme.of(context).textTheme.titleMedium,
              ),
            ),
            Padding(
              padding: const EdgeInsets.symmetric(vertical: 8.0),
              child: Text(
                'SubShape: ${shape['subShape'] ?? ''}',
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
      List<Offset> listPoints = shape['listPoints'];

      for (int i = 0; i < listPoints.length - 1; i++) {
        canvas.drawLine(listPoints[i], listPoints[i + 1], shapePaint);
      }
    } else if (shape['subShape'] == 'Horizontal Ellipse') {
      List<Offset> listPoints = shape['listPoints'];

      for (int i = 0; i < listPoints.length - 1; i++) {
        canvas.drawLine(listPoints[i], listPoints[i + 1], shapePaint);
      }
    } else if (shape['subShape'] == 'Tilted Ellipse') {
      List<Offset> listPoints = shape['listPoints'];

      for (int i = 0; i < listPoints.length - 1; i++) {
        canvas.drawLine(listPoints[i], listPoints[i + 1], shapePaint);
      }
      
    } else if (shape['shape'] == 'Triangle') {
      List<Offset> vertices = finalTriangle(shape['points']);

      for (int i = 0; i < 3; i++) {
        canvas.drawLine(vertices[i], vertices[(i + 1) % 3], shapePaint);
      }
    } else if (shape['shape'] == 'Quadrilateral' &&
        shape['subShape'] != 'Concave Quadrilateral') {
      Map<String, dynamic> verticesQuadType;
      verticesQuadType = finalQuadrilateral(shape['points']);
      List<Offset> vertices = verticesQuadType['vertices'] as List<Offset>;

      for (int i = 0; i < 4; i++) {
        canvas.drawLine(vertices[i], vertices[(i + 1) % 4], shapePaint);
      }
    } else if (shape['shape'] == 'PolyLine' ||
        shape['subShape'] == 'Concave Quadrilateral') {
      List<Offset> vertices = shape['points'];

      for (int i = 0; i < vertices.length - 1; i++) {
        canvas.drawLine(
            vertices[i], vertices[(i + 1) % vertices.length], shapePaint);
      }
    } else if (shape['shape'] == 'Convex Polygon') {
      List<Offset> polyVertices = finalPolygon(shape['points']);
      for (int i = 0; i < polyVertices.length; i++) {
        canvas.drawLine(polyVertices[i],
            polyVertices[(i + 1) % polyVertices.length], shapePaint);
      }
    } else if (shape['shape'] == 'Concave Polygon') {
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
