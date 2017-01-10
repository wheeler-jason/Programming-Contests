public class Vector2d {
	double x, y;

	public Vector2d(double x, double y) {
		this.x = x;
		this.y = y;
	}

	public Vector2d conj() {
		return new Vector2d(x, -y);
	}

	public Vector2d sub(Vector2d b) {
		return new Vector2d(x - b.x, y - b.y);
	}

	public Vector2d add(Vector2d b) {
		return new Vector2d(x + b.x, y + b.y);
	}

	public Vector2d mul(Vector2d b) {
		return new Vector2d(x * b.x - y * b.y, x * b.y + y * b.x);
	}

	public Vector2d div(Vector2d b) {
		return this.mul(b.conj()).mul(1 / b.len2());
	}

	public Vector2d mul(double b) {
		return new Vector2d(x * b, y * b);
	}

	double len2() {
		return x * x + y * y;
	}

	double len() {
		return Math.sqrt(x * x + y * y);
	}

	public Vector2d norm() {
		return len() == 0 ? new Vector2d(0, 0) : mul(1 / len());
	}

	public double cross(Vector2d b) {
		return x * b.y - y * b.x;
	}

	public double dot(Vector2d b) {
		return x * b.x + y * b.y;
	}

	public Vector2d rot() {
		return new Vector2d(-y, x);
	}

	public double proj(Vector2d p) {
		return dot(p) / len();
	}

	public static Vector2d polar(double r, double theta) {
		return new Vector2d(r * Math.cos(theta), r * Math.sin(theta));
	}

	public static Vector2d exp(Vector2d a) {
		return polar(Math.exp(a.x), a.y);
	}

	public Vector2d rotate(Vector2d p, double angle) {
		return p.sub(this).mul(exp(new Vector2d(0, angle))).add(this);
	}

	Vector2d rotate2(Vector2d p, double angle) {
		p = p.sub(this);
		double cs = Math.cos(angle);
		double sn = Math.sin(angle);
		return new Vector2d(p.x * cs - p.y * sn, p.x * sn + p.y * cs).add(this);
	}

	public Vector2d reflect(Vector2d p, Vector2d q) {
		Vector2d s = q.sub(p);
		return this.sub(p).div(s).conj().mul(s).add(p);
	}

	@Override
	public String toString() {
		return "Vector2d [x=" + x + ", y=" + y + "]";
	}

	// Usage example
	public static void main(String[] args) {
		Vector2d u = new Vector2d(0, 0);
		Vector2d v = new Vector2d(1, 0);
		Vector2d a = u.rotate(v, Math.PI * 1.0);
		Vector2d b = v.rot().rot();
		System.out.println(a);
		System.out.println(b);
	}
}

//------------------------------------------------------------------------------

public class LineGeometry {

	static final double EPS = 1e-9;

	public static int sign(double a) {
		return a < -EPS ? -1 : a > EPS ? 1 : 0;
	}

	public static class Point implements Comparable<Point> {
		public double x, y;

		public Point(double x, double y) {
			this.x = x;
			this.y = y;
		}

		public Point minus(Point b) {
			return new Point(x - b.x, y - b.y);
		}

		public double cross(Point b) {
			return x * b.y - y * b.x;
		}

		public double dot(Point b) {
			return x * b.x + y * b.y;
		}

		public Point rotateCCW(double angle) {
			return new Point(x * Math.cos(angle) - y * Math.sin(angle), x * 
                                        Math.sin(angle) + y * Math.cos(angle));
		}

		@Override
		public int compareTo(Point o) {
			// return Double.compare(Math.atan2(y, x), Math.atan2(o.y, o.x));
			return Double.compare(x, o.x) != 0 ? 
                                Double.compare(x, o.x) : Double.compare(y, o.y);
		}
	}

	public static class Line {
	    public double a, b, c;

	    public Line(double a, double b, double c) {
		    this.a = a;
		    this.b = b;
		    this.c = c;
	    }

	    public Line(Point p1, Point p2) {
		    a = +(p1.y - p2.y);
		    b = -(p1.x - p2.x);
		    c = p1.x * p2.y - p2.x * p1.y;
	    }

	    public Point intersect(Line line) {
		    double d = a * line.b - line.a * b;
		    if (sign(d) == 0) {
			    return null;
		    }
		    double x = -(c * line.b - line.c * b) / d;
		    double y = -(a * line.c - line.a * c) / d;
		    return new Point(x, y);
	    }
	}

	// Returns -1 for clockwise, 0 for straight line, 1 for ccw order
	public static int orientation(Point a, Point b, Point c) {
		Point AB = b.minus(a);
		Point AC = c.minus(a);
		return sign(AB.cross(AC));
	}

	public static boolean cw(Point a, Point b, Point c) {
		return orientation(a, b, c) < 0;
	}

	public static boolean ccw(Point a, Point b, Point c) {
		return orientation(a, b, c) > 0;
	}

	public static boolean isCrossIntersect(Point a, Point b, Point c, Point d) {
		return orientation(a, b, c) * orientation(a, b, d) < 0 &&
                                orientation(c, d, a) * orientation(c, d, b) < 0;
	}

	public static boolean isCrossOrTouchIntersect(Point a, Point b, Point c, 
                                                                    Point d) {
		if (Math.max(a.x, b.x) < Math.min(c.x, d.x) - EPS
                || Math.max(c.x, d.x) < Math.min(a.x, b.x) - EPS
				|| Math.max(a.y, b.y) < Math.min(c.y, d.y) - EPS 
                || Math.max(c.y, d.y) < Math.min(a.y, b.y) - EPS) {
			return false;
		}
		return orientation(a, b, c) * orientation(a, b, d) <= 0 
                        && orientation(c, d, a) * orientation(c, d, b) <= 0;
	}

	public static double pointToLineDistance(Point p, Line line) {
		return Math.abs(line.a * p.x + line.b * p.y + line.c) 
                                                / fastHypot(line.a, line.b);
	}

	public static double fastHypot(double x, double y) {
		return Math.sqrt(x * x + y * y);
	}

	public static double sqr(double x) {
		return x * x;
	}

	public static double angleBetween(Point a, Point b) {
		return Math.atan2(a.cross(b), a.dot(b));
	}

	public static double angle(Line line) {
		return Math.atan2(-line.a, line.b);
	}

	public static double signedArea(Point[] points) {
		int n = points.length;
		double area = 0;
		for (int i = 0, j = n - 1; i < n; j = i++) {
			area += (points[i].x - points[j].x) * (points[i].y + points[j].y);
			// area += points[i].x * points[j].y - points[j].x * points[i].y;
		}
		return area / 2;
	}

	public static enum Position {
		LEFT, RIGHT, BEHIND, BEYOND, ORIGIN, DESTINATION, BETWEEN
	}

	// Classifies position of point p against vector a
	public static Position classify(Point p, Point a) {
		int s = sign(a.cross(p));
		if (s > 0) {
			return Position.LEFT;
		}
		if (s < 0) {
			return Position.RIGHT;
		}
		if (sign(p.x) == 0 && sign(p.y) == 0) {
			return Position.ORIGIN;
		}
		if (sign(p.x - a.x) == 0 && sign(p.y - a.y) == 0) {
			return Position.DESTINATION;
		}
		if (a.x * p.x < 0 || a.y * p.y < 0) {
			return Position.BEYOND;
		}
		if (a.x * a.x + a.y * a.y < p.x * p.x + p.y * p.y) {
			return Position.BEHIND;
		}
		return Position.BETWEEN;
	}
}

//------------------------------------------------------------------------------

import java.awt.geom.*;
public class SegmentsIntersection {

	public static boolean isCrossIntersect(long x1, long y1, long x2, long y2, 
                                        long x3, long y3, long x4, long y4) {
		long z1 = (x2 - x1) * (y3 - y1) - (y2 - y1) * (x3 - x1);
		long z2 = (x2 - x1) * (y4 - y1) - (y2 - y1) * (x4 - x1);
		long z3 = (x4 - x3) * (y1 - y3) - (y4 - y3) * (x1 - x3);
		long z4 = (x4 - x3) * (y2 - y3) - (y4 - y3) * (x2 - x3);
		return (z1 < 0 || z2 < 0) && (z1 > 0 || z2 > 0) && (z3 < 0 || z4 < 0)
                                                        && (z3 > 0 || z4 > 0);
	}

	public static boolean isCrossOrTouchIntersect(long x1, long y1, long x2, 
                                long y2, long x3, long y3, long x4, long y4) {
		if (Math.max(x1, x2) < Math.min(x3, x4) 
            || Math.max(x3, x4) < Math.min(x1, x2)
		    || Math.max(y1, y2) < Math.min(y3, y4) 
            || Math.max(y3, y4) < Math.min(y1, y2)) { return false; }
	
		long z1 = (x2 - x1) * (y3 - y1) - (y2 - y1) * (x3 - x1);
		long z2 = (x2 - x1) * (y4 - y1) - (y2 - y1) * (x4 - x1);
		long z3 = (x4 - x3) * (y1 - y3) - (y4 - y3) * (x1 - x3);
		long z4 = (x4 - x3) * (y2 - y3) - (y4 - y3) * (x2 - x3);
		return (z1 <= 0 || z2 <= 0) && (z1 >= 0 || z2 >= 0) 
                && (z3 <= 0 || z4 <= 0) && (z3 >= 0 || z4 >= 0);
	}

	public static Point2D.Double getLinesIntersection(long x1, long y1, long x2,
                                 long y2, long x3, long y3, long x4, long y4) {
		long a1 = y2 - y1;
		long b1 = x1 - x2;
		long c1 = -(x1 * y2 - x2 * y1);
		long a2 = y4 - y3;
		long b2 = x3 - x4;
		long c2 = -(x3 * y4 - x4 * y3);
		long det = a1 * b2 - a2 * b1;
		if (det == 0)
			return null;
		double x = -(c1 * b2 - c2 * b1) / (double) det;
		double y = -(a1 * c2 - a2 * c1) / (double) det;
		return new Point2D.Double(x, y);
	}

	// Can also use Java API
    // Line2D.linesIntersect(x1, y1, x2, y2, x3, y3, x4, y4);
}

//------------------------------------------------------------------------------

public static double pointToSegmentDistance(int x, int y, int x1, int y1, 
    int x2, int y2) {

		long dx = x2 - x1;
		long dy = y2 - y1;
		long px = x - x1;
		long py = y - y1;
		long squaredLength = dx * dx + dy * dy;
		long dotProduct = dx * px + dy * py;
		if (dotProduct <= 0 || squaredLength == 0)
			return fastHypot(px, py);
		if (dotProduct >= squaredLength)
			return fastHypot(px - dx, py - dy);
		double q = (double) dotProduct / squaredLength;
		return fastHypot(px - q * dx, py - q * dy);
}

// Line2D.ptLineDist -- should be same result as above method within epsilon
public static double pointToLineDistance(long x, long y, long a, long b, long c) 
{
	return Math.abs(a * x + b * y + c) / fastHypot(a, b);
}

public static double fastHypot(double x, double y) {
	return Math.sqrt(x * x + y * y);
}

//------------------------------------------------------------------------------

public static Point[] convexHull(Point[] points) {
		Arrays.sort(points, (a, b) -> Integer.compare(a.x, b.x) != 0 ? 
                        Integer.compare(a.x, b.x) : Integer.compare(a.y, b.y));
		int n = points.length;
		Point[] hull = new Point[n + 1];
		int cnt = 0;
		for (int i = 0; i < 2 * n; i++) {
			int j = i < n ? i : 2 * n - 1 - i;
			while (cnt >= 2 && 
                       isNotRightTurn(hull[cnt - 2], hull[cnt - 1], points[j]))
				--cnt;

			hull[cnt++] = points[j];
		}
		return Arrays.copyOf(hull, cnt - 1);
}

static boolean isNotRightTurn(Point a, Point b, Point c) {
	long cross = (long)(a.x - b.x)*(c.y - b.y) - (long)(a.y - b.y)*(c.x - b.x);
	long dot = (long)(a.x - b.x)*(c.x - b.x) + (long)(a.y - b.y)*(c.y - b.y);
	return cross < 0 || cross == 0 && dot <= 0;
}

public static class Point {
	public final int x, y;

	public Point(int x, int y) {
		this.x = x;
		this.y = y;
	}
}
//------------------------------------------------------------------------------

// qx and qy denote the point 
// returns -1 if outside, 0 if on boundary, 1 if inside
public static int pointInPolygon(int qx, int qy, int[] x, int[] y) {
	int n = x.length;
	int cnt = 0;
	for (int i = 0, j = n - 1; i < n; j = i++) {
		if (y[i] == qy && (x[i] == qx || y[j] == qy && (x[i] <= qx 
            || x[j] <= qx) && (x[i] >= qx || x[j] >= qx))) { return 0; }

		if ((y[i] > qy) != (y[j] > qy)) {
			long det = ((long) x[i] - qx) * ((long) y[j] - qy) 
                        - ((long) x[j] - qx) * ((long) y[i] - qy);
			if (det == 0)
				return 0; 
			if ((det > 0) != (y[j] > y[i]))
				++cnt;
		}
	}
	return cnt % 2 == 0 ? -1 : 1;
}

//------------------------------------------------------------------------------

public class CircleOperations {

	static final double EPS = 1e-9;

	public static double fastHypot(double x, double y) {
		return Math.sqrt(x * x + y * y);
	}

	public static class Point {
		double x, y;

		public Point(double x, double y) {
			this.x = x;
			this.y = y;
		}
	}

	public static class Circle {
		double x, y, r;

		public Circle(double x, double y, double r) {
			this.x = x;
			this.y = y;
			this.r = r;
		}

		public boolean contains(Point p) {
			return fastHypot(p.x - x, p.y - y) < r + EPS;
		}
	}

	public static class Line {
		double a, b, c;

		public Line(double a, double b, double c) {
			this.a = a;
			this.b = b;
			this.c = c;
		}

		public Line(Point p1, Point p2) {
			a = +(p1.y - p2.y);
			b = -(p1.x - p2.x);
			c = p1.x * p2.y - p2.x * p1.y;
		}
	}

	public static Point[] circleLineIntersection(Circle circle, Line line) {
		return Math.abs(line.a) >= Math.abs(line.b)
		? intersection(line.a,line.b,line.c,circle.x,circle.y,circle.r,false)
		: intersection(line.b,line.a,line.c,circle.y,circle.x,circle.r,true);
	}

	static Point[] intersection(double a, double b, double c, double CX, 
                                          double CY, double R, boolean swap) {
		double A = a * a + b * b;
		double B = 2.0 * b * (c + a * CX) - 2.0 * a * a * CY;
		double C = (c + a * CX) * (c + a * CX) + a * a * (CY * CY - R * R);
		double d = B * B - 4 * A * C;
		if (d < -EPS)
			return new Point[0];
		d = Math.sqrt(d < 0 ? 0 : d);
		double y1 = (-B + d) / (2 * A);
		double x1 = (-c - b * y1) / a;
		double y2 = (-B - d) / (2 * A);
		double x2 = (-c - b * y2) / a;
		return swap ? d > EPS ? new Point[]{new Point(y1, x1), new Point(y2, x2)}
                              : new Point[]{new Point(y1, x1)}
				    : d > EPS ? new Point[]{new Point(x1, y1), new Point(x2, y2)}
                              : new Point[]{new Point(x1, y1)};
	}

	public static Point[] circleCircleIntersection(Circle c1, Circle c2) {
		if (fastHypot(c1.x - c2.x, c1.y - c2.y) < EPS) {
			if (Math.abs(c1.r - c2.r) < EPS)
				return null; // infinity intersections
			return new Point[0];
		}
		double dx = c2.x - c1.x;
		double dy = c2.y - c1.y;
		double A = -2 * dx;
		double B = -2 * dy;
		double C = dx * dx + dy * dy + c1.r * c1.r - c2.r * c2.r;
		Point[] res = circleLineIntersection(new Circle(0, 0, c1.r), 
                                                             new Line(A, B, C));
		for (Point point : res) {
			point.x += c1.x;
			point.y += c1.y;
		}
		return res;
	}

	public static double circleCircleIntersectionArea(Circle c1, Circle c2) {
		double r = Math.min(c1.r, c2.r);
		double R = Math.max(c1.r, c2.r);
		double d = fastHypot(c1.x - c2.x, c1.y - c2.y);
		if (d < R - r + EPS)
			return Math.PI * r * r;
		if (d > R + r - EPS)
			return 0;
		double area = r*r*Math.acos((d * d + r * r - R * R) / 2 / d / r) + R * R
				* Math.acos((d * d + R * R - r * r) / 2 / d / R) - 0.5
				* Math.sqrt((-d + r + R)*(d + r - R)*(d - r + R)*(d + r + R));
		return area;
	}

	public static Circle minEnclosingCircle(Point[] points) {
		if (points.length == 0)
			return new Circle(0, 0, 0);
		if (points.length == 1)
			return new Circle(points[0].x, points[0].y, 0);
		Collections.shuffle(Arrays.asList(points));
		Circle circle = getCircumCircle(points[0], points[1]);
		for (int i = 2; i < points.length; i++) {
			if (!circle.contains(points[i])) {
				circle = getCircumCircle(points[0], points[i]);
				for (int j = 1; j < i; j++) {
					if (!circle.contains(points[j])) {
						circle = getCircumCircle(points[j], points[i]);
						for (int k = 0; k < j; k++) {
							if (!circle.contains(points[k])) {
								circle = 
                                 getCircumCircle(points[i],points[j],points[k]);
							}
						}
					}
				}
			}
		}
		return circle;
	}

	public static Circle getCircumCircle(Point a, Point b) {
		double x = (a.x + b.x) / 2.;
		double y = (a.y + b.y) / 2.;
		double r = fastHypot(a.x - x, a.y - y);
		return new Circle(x, y, r);
	}

	public static Circle getCircumCircle(Point a, Point b, Point c) {
		double Bx = b.x - a.x;
		double By = b.y - a.y;
		double Cx = c.x - a.x;
		double Cy = c.y - a.y;
		double d = 2 * (Bx * Cy - By * Cx);
		if (Math.abs(d) < EPS)
			return getCircumCircle(new Point(Math.min(a.x, Math.min(b.x, c.x)), 
                    Math.min(a.y, Math.min(b.y, c.y))),
					new Point(Math.max(a.x, Math.max(b.x, c.x)), 
                    Math.max(a.y, Math.max(b.y, c.y))));
		double z1 = Bx * Bx + By * By;
		double z2 = Cx * Cx + Cy * Cy;
		double cx = Cy * z1 - By * z2;
		double cy = Bx * z2 - Cx * z1;
		double x = cx / d;
		double y = cy / d;
		double r = fastHypot(x, y);
		return new Circle(x + a.x, y + a.y, r);
	}

	static boolean eq(Point p1, Point p2) {
		return !(fastHypot(p1.x - p2.x, p1.y - p2.y) > EPS);
	}
}

//------------------------------------------------------------------------------

