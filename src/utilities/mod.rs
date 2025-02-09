pub fn rescale<T: Copy + std::ops::Add<Output = T> + std::ops::Sub<Output = T> + std::ops::Div<Output = T>+ std::ops::Mul<Output = T>>
(value: T, old_min: T, old_max: T, new_min: T, new_max: T) -> T {
    (((value - old_min) / (old_max - old_min)) * (new_max - new_min)) + new_min
}

// copied from fft_2d because it is not a public function:
pub fn transpose<T: Copy + Default>(width: usize, height: usize, matrix: &[T]) -> Vec<T> {
    let mut ind = 0;
    let mut ind_tr;
    let mut transposed = vec![T::default(); matrix.len()];
    for row in 0..height {
        ind_tr = row;
        for _ in 0..width {
            transposed[ind_tr] = matrix[ind];
            ind += 1;
            ind_tr += height;
        }
    }
    transposed
}