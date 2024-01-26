use std::{
    fmt::Display,
    ops::{Index, IndexMut},
};

#[derive(Debug, PartialEq)]
pub struct Matrix<T> {
    rows: usize,
    cols: usize,
    data: Box<[T]>,
}

impl<T: Default + Clone> Matrix<T> {
    pub fn zeros(rows: usize, cols: usize) -> Self {
        Self {
            rows,
            cols,
            data: vec![T::default(); rows * cols].into_boxed_slice(),
        }
    }
}

impl<T> Matrix<T> {
    pub fn new(v: Vec<Vec<T>>) -> Self {
        let rows = v.len();
        let v: Vec<T> = v.into_iter().flatten().collect();
        let cols = v.len() / rows;
        Self {
            rows,
            cols,
            data: v.into_boxed_slice(),
        }
    }

    pub fn data(&self) -> &[T] {
        &self.data
    }

    /// (rows, cols)
    pub fn shape(&self) -> (usize, usize) {
        (self.rows, self.cols)
    }
}

impl<T> Index<(usize, usize)> for Matrix<T> {
    type Output = T;

    fn index(&self, (x, y): (usize, usize)) -> &Self::Output {
        &self.data[self.cols * x + y]
    }
}

impl<T> IndexMut<(usize, usize)> for Matrix<T> {
    fn index_mut(&mut self, (x, y): (usize, usize)) -> &mut Self::Output {
        &mut self.data[self.cols * x + y]
    }
}

impl<T: Display> Display for Matrix<T> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let width = f.width().unwrap_or(8);
        let prec = f.precision().unwrap_or(4);
        for row in 0..self.rows {
            for col in 0..self.cols {
                write!(f, "{col:width$.prec$}", col = &self[(row, col)])?;
            }
            writeln!(f)?;
        }
        Ok(())
    }
}
