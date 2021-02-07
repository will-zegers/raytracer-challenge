use std::cmp::PartialEq;
use std::f64::consts::PI;
use std::ops::{Index, IndexMut, Mul};

use crate::point::Point;
use crate::vector::Vector;

const TOL: f64 = 1e-6;

pub enum Axis {
    X,
    Y,
    Z,
}

#[derive(Clone, Debug)]
pub struct Matrix {
    m: usize,
    n: usize,
    elements: Vec<f64>,
}

impl Matrix {
    pub fn new(m: usize, n: usize, elements: Vec<f64>) -> Self {
        assert_eq!(elements.len(), m * n);
        Self { m, n, elements }
    }

    pub fn eye(n: usize) -> Self {
        let mut elements = vec![0.; n * n];
        for i in 0..n {
            for j in 0..n {
                if i == j {
                    elements[i * n + j] = 1.
                }
            }
        }
        Self { m: n, n, elements }
    }

    pub fn transpose(mut self) -> Self {
        for i in 0..self.m {
            for j in i..self.n {
                let idx1 = i * self.n + j;
                let idx2 = j * self.n + i;
                self.elements.swap(idx1, idx2);
            }
        }
        self
    }

    pub fn det2x2(&self) -> f64 {
        debug_assert!(self.m == 2);
        debug_assert!(self.n == 2);

        self.elements[0] * self.elements[3] - self.elements[1] * self.elements[2]
    }

    pub fn submat(&self, r: usize, c: usize) -> Self {
        let m = self.m - 1;
        let n = self.n - 1;
        let mut v = Vec::with_capacity(m * n);

        for i in 0..self.m {
            if i == r {
                continue;
            }
            for j in 0..self.n {
                if j != c {
                    v.push(self.elements[i * self.n + j]);
                }
            }
        }
        Self { m, n, elements: v }
    }

    pub fn minor(&self, r: usize, c: usize) -> f64 {
        self.submat(r, c).det()
    }

    pub fn cofactor(&self, r: usize, c: usize) -> f64 {
        match (r + c) % 2 {
            0 => self.minor(r, c),
            _ => -self.minor(r, c),
        }
    }

    pub fn det(&self) -> f64 {
        if self.m != self.n {
            return 0.;
        }

        if self.n == 2 {
            return self.det2x2();
        }

        let mut d = 0.;
        for i in 0..self.n {
            d += self.elements[i] * self.cofactor(0, i);
        }

        d
    }

    fn is_invertible(&self) -> bool {
        return self.det() != 0.;
    }

    pub fn inverse(&self) -> Self {
        assert!(self.is_invertible());

        let mut v = Vec::with_capacity(self.m * self.n);
        let d = self.det();
        for i in 0..self.m {
            for j in 0..self.n {
                v.push(self.cofactor(i, j) / d);
            }
        }
        Self::new(self.m, self.n, v).transpose()
    }
}

// transform factories
impl Matrix {
    pub fn translation(x: f64, y: f64, z: f64) -> Self {
        Self {
            m: 4,
            n: 4,
            elements: vec![1., 0., 0., x, 0., 1., 0., y, 0., 0., 1., z, 0., 0., 0., 1.],
        }
    }

    pub fn scaling(x: f64, y: f64, z: f64) -> Self {
        Self {
            m: 4,
            n: 4,
            elements: vec![x, 0., 0., 0., 0., y, 0., 0., 0., 0., z, 0., 0., 0., 0., 1.],
        }
    }

    pub fn rotation(axis: Axis, r: f64) -> Self {
        let mut m = Matrix::eye(4);
        match axis {
            Axis::X => {
                m[(1, 1)] = f64::cos(r);
                m[(1, 2)] = -f64::sin(r);
                m[(2, 1)] = f64::sin(r);
                m[(2, 2)] = f64::cos(r);
            }
            Axis::Y => {
                m[(0, 0)] = f64::cos(r);
                m[(0, 2)] = f64::sin(r);
                m[(2, 0)] = -f64::sin(r);
                m[(2, 2)] = f64::cos(r);
            }
            Axis::Z => {
                m[(0, 0)] = f64::cos(r);
                m[(0, 1)] = -f64::sin(r);
                m[(1, 0)] = f64::sin(r);
                m[(1, 1)] = f64::cos(r);
            }
        };
        m
    }

    fn shearing(xy: f64, xz: f64, yx: f64, yz: f64, zx: f64, zy: f64) -> Self {
        Self {
            m: 4,
            n: 4,
            elements: vec![
                1., xy, xz, 0., yx, 1., yz, 0., zx, zy, 1., 0., 0., 0., 0., 1.,
            ],
        }
    }
}

impl PartialEq<Matrix> for Matrix {
    fn eq(&self, rhs: &Matrix) -> bool {
        for (e1, e2) in self.elements.iter().zip(rhs.elements.iter()) {
            if f64::abs(e1 - e2) > TOL {
                return false;
            }
        }
        true
    }
}

impl Index<(usize, usize)> for Matrix {
    type Output = f64;

    fn index(&self, i: (usize, usize)) -> &Self::Output {
        &self.elements[i.0 * self.m + i.1]
    }
}

impl IndexMut<(usize, usize)> for Matrix {
    fn index_mut(&mut self, i: (usize, usize)) -> &mut Self::Output {
        &mut self.elements[i.0 * self.m + i.1]
    }
}

impl Mul<Matrix> for Matrix {
    type Output = Matrix;

    fn mul(self, rhs: Matrix) -> Self::Output {
        assert_eq!(self.n, rhs.m);

        let mut v = vec![0.; self.m * rhs.n];
        for i in 0..self.m {
            for j in 0..rhs.n {
                for k in 0..self.n {
                    let idx1 = i * self.n + k;
                    let idx2 = k * rhs.n + j;
                    let idx3 = i * rhs.n + j;
                    v[idx3] += self.elements[idx1] * rhs.elements[idx2];
                }
            }
        }

        Self::Output {
            m: self.m,
            n: rhs.n,
            elements: v,
        }
    }
}

impl Mul<Point> for &Matrix {
    type Output = Point;

    fn mul(self, rhs: Point) -> Self::Output {
        let mut v = vec![0.; self.m];
        for i in 0..self.m {
            for j in 0..self.n {
                let idx1 = i * self.n + j;
                v[i] += self.elements[idx1] * rhs[j];
            }
        }

        Self::Output::new(v[0], v[1], v[2])
    }
}

impl Mul<Point> for Matrix {
    type Output = Point;

    fn mul(self, rhs: Point) -> Self::Output {
        &self * rhs
    }
}

impl Mul<Vector> for &Matrix {
    type Output = Vector;

    fn mul(self, rhs: Vector) -> Self::Output {
        let mut v = vec![0.; self.m];
        for i in 0..self.m {
            for j in 0..self.n {
                let idx1 = i * self.n + j;
                v[i] += self.elements[idx1] * rhs[j];
            }
        }

        Self::Output::new(v[0], v[1], v[2])
    }
}

impl Mul<Vector> for Matrix {
    type Output = Vector;

    fn mul(self, rhs: Vector) -> Self::Output {
        &self * rhs
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn new() {
        let v = vec![-3., 5., 1., -2.];
        let m = Matrix::new(2, 2, v);

        assert_eq!(m[(0, 0)], -3.);
        assert_eq!(m[(0, 1)], 5.);
        assert_eq!(m[(1, 0)], 1.);
        assert_eq!(m[(1, 1)], -2.);

        let v = vec![-3., 5., 0., 1., -2., -7., 0., 1., 1.];
        let m = Matrix::new(3, 3, v);
        assert_eq!(m[(0, 0)], -3.);
        assert_eq!(m[(1, 1)], -2.);
        assert_eq!(m[(2, 2)], 1.);

        let v = vec![
            1., 2., 3., 4., 5.5, 6.5, 7.5, 8.5, 9., 10., 11., 12., 13.5, 14.5, 15.5, 16.5,
        ];
        let m = Matrix::new(4, 4, v);

        assert_eq!(m[(0, 0)], 1.);
        assert_eq!(m[(0, 3)], 4.);
        assert_eq!(m[(1, 0)], 5.5);
        assert_eq!(m[(1, 2)], 7.5);
        assert_eq!(m[(2, 2)], 11.);
        assert_eq!(m[(3, 0)], 13.5);
        assert_eq!(m[(3, 2)], 15.5);
    }

    #[test]
    fn eq() {
        let v1 = vec![
            1., 2., 3., 4., 5., 6., 7., 8., 9., 8., 7., 6., 5., 4., 3., 2.,
        ];
        let m1 = Matrix::new(4, 4, v1);

        let v2 = vec![
            1., 2., 3., 4., 5., 6., 7., 8., 9., 8., 7., 6., 5., 4., 3., 2.,
        ];
        let m2 = Matrix::new(4, 4, v2);

        assert_eq!(m1, m2);

        let v3 = vec![
            2., 3., 4., 5., 6., 7., 8., 9., 8., 7., 6., 5., 4., 3., 2., 1.,
        ];
        let m3 = Matrix::new(4, 4, v3);

        assert_ne!(m1, m3);
    }

    #[test]
    fn mul() {
        // Matrix-matrix
        let v1 = vec![
            1., 2., 3., 4., 5., 6., 7., 8., 9., 8., 7., 6., 5., 4., 3., 2.,
        ];
        let m1 = Matrix::new(4, 4, v1);

        let v2 = vec![
            -2., 1., 2., 3., 3., 2., 1., -1., 4., 3., 6., 5., 1., 2., 7., 8.,
        ];
        let m2 = Matrix::new(4, 4, v2);

        assert_eq!(
            m1 * m2,
            Matrix::new(
                4,
                4,
                vec![
                    20.0, 22.0, 50.0, 48.0, 44.0, 54.0, 114.0, 108.0, 40.0, 58.0, 110.0, 102.0,
                    16.0, 26.0, 46.0, 42.0
                ]
            )
        );

        // Matrix-point
        let a = Matrix::new(
            4,
            4,
            vec![
                1., 2., 3., 4., 2., 4., 4., 2., 8., 6., 4., 1., 0., 0., 0., 1.,
            ],
        );
        let b = Point::new(1., 2., 3.);
        assert_eq!(a * b, Point::new(18., 24., 33.));

        // Matrix-vector
        let a = Matrix::new(
            4,
            4,
            vec![
                1., 2., 3., 4., 2., 4., 4., 2., 8., 6., 4., 1., 0., 0., 0., 1.,
            ],
        );
        let b = Vector::new(1., 2., 3.);
        assert_eq!(a * b, Vector::new(14., 22., 32.));
    }

    #[test]
    fn eye() {
        let m = Matrix::new(
            4,
            4,
            vec![
                0., 1., 2., 4., 1., 2., 4., 8., 2., 4., 8., 16., 4., 8., 16., 32.,
            ],
        );
        assert_eq!(m.clone() * Matrix::eye(4), m);

        let p = Point::new(1., 2., 3.);
        assert_eq!(Matrix::eye(3) * p.clone(), p);
    }

    #[test]
    fn transpose() {
        let m = Matrix::new(
            4,
            4,
            vec![
                0., 9., 3., 0., 9., 8., 0., 8., 1., 8., 5., 3., 0., 0., 5., 8.,
            ],
        );
        let m_t = Matrix::new(
            4,
            4,
            vec![
                0., 9., 1., 0., 9., 8., 8., 0., 3., 0., 5., 5., 0., 8., 3., 8.,
            ],
        );
        assert_eq!(m.transpose(), m_t);
    }

    #[test]
    fn det2x2() {
        let m = Matrix::new(2, 2, vec![1., 5., -3., 2.]);
        assert_eq!(m.det2x2(), 17.);
    }

    #[test]
    fn submat() {
        let m = Matrix::new(3, 3, vec![1., 5., 0., -3., 2., 7., 0., 6., -3.]);
        assert_eq!(m.submat(0, 2), Matrix::new(2, 2, vec![-3., 2., 0., 6.]));

        let m = Matrix::new(
            4,
            4,
            vec![
                -6., 1., 1., 6., -8., 5., 8., 6., -1., 0., 9., 2., -7., 1., -1., 1.,
            ],
        );
        assert_eq!(
            m.submat(2, 1),
            Matrix::new(3, 3, vec![-6., 1., 6., -8., 8., 6., -7., -1., 1.])
        );
    }

    #[test]
    fn minor() {
        let m = Matrix::new(3, 3, vec![3., 5., 0., 2., -1., -7., 6., -1., 5.]);
        assert_eq!(m.minor(1, 0), 25.)
    }

    #[test]
    fn cofactor() {
        let m = Matrix::new(3, 3, vec![3., 5., 0., 2., -1., -7., 6., -1., 5.]);
        assert_eq!(m.minor(0, 0), -12.);
        assert_eq!(m.cofactor(0, 0), -12.);
        assert_eq!(m.minor(1, 0), 25.);
        assert_eq!(m.cofactor(1, 0), -25.);
    }

    #[test]
    fn det() {
        let m = Matrix::new(3, 3, vec![1., 2., 6., -5., 8., -4., 2., 6., 4.]);

        assert_eq!(m.cofactor(0, 0), 56.);
        assert_eq!(m.cofactor(0, 1), 12.);
        assert_eq!(m.cofactor(0, 2), -46.);
        assert_eq!(m.det(), -196.);

        let m = Matrix::new(
            4,
            4,
            vec![
                -2., -8., 3., 5., -3., 1., 7., 3., 1., 2., -9., 6., -6., 7., 7., -9.,
            ],
        );
        assert_eq!(m.cofactor(0, 0), 690.);
        assert_eq!(m.cofactor(0, 1), 447.);
        assert_eq!(m.cofactor(0, 2), 210.);
        assert_eq!(m.cofactor(0, 3), 51.);
        assert_eq!(m.det(), -4071.);
    }

    #[test]
    fn is_invertible() {
        let m = Matrix::new(
            4,
            4,
            vec![
                6., 4., 4., 4., 5., 5., 7., 6., 4., -9., 3., -7., 9., 1., 7., -6.,
            ],
        );
        assert_eq!(m.det(), -2120.);
        assert!(m.is_invertible());

        let m = Matrix::new(
            4,
            4,
            vec![
                -4., 2., -2., -3., 9., 6., 2., 6., 0., -5., 1., -5., 0., 0., 0., 0.,
            ],
        );
        assert_eq!(m.det(), 0.);
        assert!(!m.is_invertible());
    }

    #[test]
    fn invert() {
        let m = Matrix::new(
            4,
            4,
            vec![
                -5., 2., 6., -8., 1., -5., 1., 8., 7., 7., -6., -7., 1., -3., 7., 4.,
            ],
        );
        let m_inv = m.inverse();

        assert_eq!(m.det(), 532.);
        assert_eq!(m.cofactor(2, 3), -160.);
        assert_eq!(m_inv[(3, 2)], -160. / 532.);

        assert_eq!(m.cofactor(3, 2), 105.);
        assert_eq!(m_inv[(2, 3)], 105. / 532.);

        assert_eq!(
            m_inv,
            Matrix::new(
                4,
                4,
                vec![
                    0.218045112,
                    0.451127819,
                    0.240601503,
                    -0.045112781,
                    -0.808270676,
                    -1.456766917,
                    -0.443609022,
                    0.520676691,
                    -0.078947368,
                    -0.223684210,
                    -0.052631578,
                    0.197368421,
                    -0.522556390,
                    -0.813909774,
                    -0.300751879,
                    0.306390977,
                ]
            )
        );

        let m = Matrix::new(
            4,
            4,
            vec![
                8., -5., 9., 2., 7., 5., 6., 1., -6., 0., 9., 6., -3., 0., -9., -4.,
            ],
        );
        assert_eq!(
            m.inverse(),
            Matrix::new(
                4,
                4,
                vec![
                    -0.153846, -0.153846, -0.282051, -0.538461, -0.076923, 0.123076, 0.025641,
                    0.030769, 0.358974, 0.358974, 0.435897, 0.923076, -0.692307, -0.692307,
                    -0.769230, -1.923076,
                ]
            )
        );

        let m = Matrix::new(
            4,
            4,
            vec![
                8., -5., 9., 2., 7., 5., 6., 1., -6., 0., 9., 6., -3., 0., -9., -4.,
            ],
        );
        assert_eq!(
            m.inverse(),
            Matrix::new(
                4,
                4,
                vec![
                    -0.153846, -0.153846, -0.282051, -0.538461, -0.076923, 0.1230769, 0.0256410,
                    0.0307692, 0.3589743, 0.3589743, 0.4358974, 0.9230769, -0.692307, -0.692307,
                    -0.769230, -1.923076,
                ]
            )
        );

        let m = Matrix::new(
            4,
            4,
            vec![
                9., 3., 0., 9., -5., -2., -6., -3., -4., 9., 6., 4., -7., 6., 6., 2.,
            ],
        );
        assert_eq!(
            m.inverse(),
            Matrix::new(
                4,
                4,
                vec![
                    -0.040740, -0.077777, 0.1444444, -0.222222, -0.077777, 0.0333333, 0.3666666,
                    -0.333333, -0.029012, -0.146296, -0.109259, 0.1296296, 0.1777777, 0.0666666,
                    -0.266666, 0.3333333,
                ]
            )
        );

        let m1 = Matrix::new(
            4,
            4,
            vec![
                3., -9., 7., 3., 3., -8., 2., -9., -4., 4., 4., 1., -6., 5., -1., 1.,
            ],
        );
        let m2 = Matrix::new(
            4,
            4,
            vec![
                8., 2., 2., 2., 3., -1., 7., 0., 7., 0., 5., 4., 6., -2., 0., 5.,
            ],
        );
        let m3 = m1.clone() * m2.clone();
        assert_eq!(m3 * m2.inverse(), m1);
    }

    #[test]
    fn translation() {
        // point translation
        let t = Matrix::translation(5., -3., 2.);
        let p = Point::new(-3., 4., 5.);
        assert_eq!(t * p, Point::new(2., 1., 7.));

        let t = Matrix::translation(5., -3., 2.);
        let p = Point::new(0., 0., 0.);
        assert_eq!(t * p, Point::new(5., -3., 2.));

        // vector translation
        let t = Matrix::translation(5., -3., 2.);
        let p = Vector::new(-3., 4., 5.);
        assert_eq!(t * p, p);
    }

    #[test]
    fn scaling() {
        // point scaling
        let t = Matrix::scaling(2., 3., 4.);
        let p = Point::new(-4., 6., 8.);
        assert_eq!(t * p, Point::new(-8., 18., 32.));

        // vector scaling
        let t = Matrix::scaling(2., 3., 4.);
        let v = Vector::new(-4., 6., 8.);
        assert_eq!(t * v, Vector::new(-8., 18., 32.));

        // inverse scaling
        let t = Matrix::scaling(2., 3., 4.);
        let t_inv = t.inverse();
        let v = Vector::new(-4., 6., 8.);
        assert_eq!(t_inv * v, Vector::new(-2., 2., 2.));

        // reflection (negative scaling)
        let t = Matrix::scaling(-1., 1., 1.);
        let p = Point::new(2., 3., 4.);
        assert_eq!(t * p, Point::new(-2., 3., 4.));
    }

    #[test]
    fn rotation() {
        // rotate around the x-axis
        let p = Point::new(0., 1., 0.);
        let half_quarter = Matrix::rotation(Axis::X, PI / 4.);
        let full_quarter = Matrix::rotation(Axis::X, PI / 2.);

        assert_eq!(
            half_quarter * p,
            Point::new(0., f64::sqrt(2.) / 2., f64::sqrt(2.) / 2.)
        );
        assert_eq!(full_quarter * p, Point::new(0., 0., 1.));

        // rotate around the x-axis in the opposite direction
        let p = Point::new(0., 1., 0.);
        let half_quarter = Matrix::rotation(Axis::X, PI / 4.);
        let inv = half_quarter.inverse();
        assert_eq!(
            inv * p,
            Point::new(0., f64::sqrt(2.) / 2., -f64::sqrt(2.) / 2.)
        );

        // rotate around the y-axis
        let p = Point::new(0., 0., 1.);
        let half_quarter = Matrix::rotation(Axis::Y, PI / 4.);
        let full_quarter = Matrix::rotation(Axis::Y, PI / 2.);

        assert_eq!(
            half_quarter * p,
            Point::new(f64::sqrt(2.) / 2., 0., f64::sqrt(2.) / 2.)
        );
        assert_eq!(full_quarter * p, Point::new(1., 0., 0.));

        // rotate around the z-axis
        let p = Point::new(0., 1., 0.);
        let half_quarter = Matrix::rotation(Axis::Z, PI / 4.);
        let full_quarter = Matrix::rotation(Axis::Z, PI / 2.);

        assert_eq!(
            half_quarter * p,
            Point::new(-f64::sqrt(2.) / 2., f64::sqrt(2.) / 2., 0.)
        );
        assert_eq!(full_quarter * p, Point::new(-1., 0., 0.));
    }

    #[test]
    fn shearing() {
        // shear in x in propotion to y
        let t = Matrix::shearing(1., 0., 0., 0., 0., 0.);
        let p = Point::new(2., 3., 4.);
        assert_eq!(t * p, Point::new(5., 3., 4.));

        // shear in x in propotion to z
        let t = Matrix::shearing(0., 1., 0., 0., 0., 0.);
        let p = Point::new(2., 3., 4.);
        assert_eq!(t * p, Point::new(6., 3., 4.));

        // shear in y in propotion to x
        let t = Matrix::shearing(0., 0., 1., 0., 0., 0.);
        let p = Point::new(2., 3., 4.);
        assert_eq!(t * p, Point::new(2., 5., 4.));

        // shear in y in propotion to z
        let t = Matrix::shearing(0., 0., 0., 1., 0., 0.);
        let p = Point::new(2., 3., 4.);
        assert_eq!(t * p, Point::new(2., 7., 4.));

        // shear in z in propotion to x
        let t = Matrix::shearing(0., 0., 0., 0., 1., 0.);
        let p = Point::new(2., 3., 4.);
        assert_eq!(t * p, Point::new(2., 3., 6.));

        // shear in z in propotion to y
        let t = Matrix::shearing(0., 0., 0., 0., 0., 1.);
        let p = Point::new(2., 3., 4.);
        assert_eq!(t * p, Point::new(2., 3., 7.));
    }

    #[test]
    fn chained() {
        let p0 = Point::new(1., 0., 1.);
        let m_r = Matrix::rotation(Axis::X, PI / 2.);
        let m_s = Matrix::scaling(5., 5., 5.);
        let m_t = Matrix::translation(10., 5., 7.);

        // rotate first
        let p1 = m_r * p0;
        assert_eq!(p1, Point::new(1., -1., 0.));

        // then apply scaling
        let p2 = m_s * p1;
        assert_eq!(p2, Point::new(5., -5., 0.));

        // then apply translation
        let p4 = m_t * p2;
        assert_eq!(p4, Point::new(15., 0., 7.));

        // chain multiple matrices
        let p = Point::new(1., 0., 1.);
        let m_r = Matrix::rotation(Axis::X, PI / 2.);
        let m_s = Matrix::scaling(5., 5., 5.);
        let m_t = Matrix::translation(10., 5., 7.);
        let t = m_t * m_s * m_r;

        assert_eq!(t * p, Point::new(15., 0., 7.));
    }
}
