use std::{
    fmt::Display,
    ops::{Index, IndexMut},
};

#[derive(Debug, PartialEq)]
pub struct Matrix<T>(Vec<Vec<T>>);

impl<T: Default + Clone> Matrix<T> {
    pub fn zeros(rows: usize, cols: usize) -> Self {
        Self(vec![vec![T::default(); cols]; rows])
    }
}

impl<T> Matrix<T> {
    pub fn new(v: Vec<Vec<T>>) -> Self {
        Self(v)
    }

    pub fn data(&self) -> &Vec<Vec<T>> {
        &self.0
    }

    /// (rows, cols)
    pub fn shape(&self) -> (usize, usize) {
        (self.0.len(), self.0.first().map(|v| v.len()).unwrap_or(0))
    }
}

impl<T> Index<(usize, usize)> for Matrix<T> {
    type Output = T;

    fn index(&self, (x, y): (usize, usize)) -> &Self::Output {
        &self.0[x][y]
    }
}

impl<T> IndexMut<(usize, usize)> for Matrix<T> {
    fn index_mut(&mut self, (x, y): (usize, usize)) -> &mut Self::Output {
        &mut self.0[x][y]
    }
}

impl<T: Display> Display for Matrix<T> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let width = f.width().unwrap_or(8);
        let prec = f.precision().unwrap_or(4);
        for row in &self.0 {
            for col in row {
                write!(f, "{col:width$.prec$}")?;
            }
            writeln!(f)?;
        }
        Ok(())
    }
}
