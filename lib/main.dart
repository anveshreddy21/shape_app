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
  // A list of points representing the current stroke
  List<Offset> currentStroke = [];
  List<Offset> rdpPoints = [];
  List<int> rdpIndices = [];
  List<double> localAngles = [];
  Map<String, dynamic> shape = {};

  @override
  Widget build(BuildContext context) {
    return Scaffold(
      appBar: AppBar(
        title: const Text('Drawing Shapes App'),
        backgroundColor: Color.fromARGB(255, 52, 155, 158),
      ),
      body: Builder(
        builder: (context) => GestureDetector(
          onPanUpdate: (details) {
            setState(() {
              RenderBox renderBox = context.findRenderObject() as RenderBox;
              Offset localPosition =
                  renderBox.globalToLocal(details.globalPosition);
              currentStroke.add(localPosition);
            });
          },
          onPanDown: (details) {
            // Clear the current stroke and start a new one
            setState(() {
              currentStroke.clear();
              rdpPoints.clear();
              rdpIndices.clear();
              localAngles.clear();
              shape.clear();
              RenderBox renderBox = context.findRenderObject() as RenderBox;
              Offset localPosition =
                  renderBox.globalToLocal(details.globalPosition);
              currentStroke.add(localPosition);
            });
          },
          onPanEnd: (details) {
            currentStroke = closure(currentStroke);
            var result = rdpNRWithIndices(currentStroke, alpha: 0.05);

            setState(() {
              rdpPoints = result['points'] as List<Offset>;
              rdpIndices = result['indices'] as List<int>;
              localAngles = calculateLocalAngles(currentStroke, rdpIndices);
              shape = shapeDecider(currentStroke, rdpPoints, localAngles);
             
            });
          },
          child: CustomPaint(
            painter: MyPainter(currentStroke, rdpPoints, shape),
            size: Size.infinite,
          ),
        ),
      ),
    );
  }
}

class MyPainter extends CustomPainter {
  final List<Offset> points;
  final List<Offset> rdpPoints;
  final Map<String, dynamic> shape;

  MyPainter(this.points, this.rdpPoints, this.shape);

  @override
  void paint(Canvas canvas, Size size) {
    Paint paint = Paint()
      ..color = Colors.black
      ..strokeCap = StrokeCap.round
      ..strokeWidth = 6.0;

    for (int i = 0; i < points.length - 1; i++) {
      canvas.drawLine(points[i], points[i + 1], paint);
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


