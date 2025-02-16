use num_traits::{Num, pow};
use rustfft::FftPlanner;
use rustfft::num_complex::Complex;use std::fs;
use std::path::Path;

/// Ensures the specified directory exists, creating it if necessary.
pub fn ensure_directory_exists(dir: &str) -> std::io::Result<()> {
    let path = Path::new(dir);
    fs::create_dir_all(path)?;
    Ok(())
}

/// Rescales a value from one range to another.
///
/// This function maps a `value` from the range `[old_min, old_max]` to a new range `[new_min, new_max]`.
///
/// # Parameters:
/// - `value`: The value to be rescaled.
/// - `old_min`: The minimum of the original range.
/// - `old_max`: The maximum of the original range.
/// - `new_min`: The minimum of the target range.
/// - `new_max`: The maximum of the target range.
///
/// # Returns:
/// The rescaled value within the new range.
pub fn rescale<T>(value: T, old_min: T, old_max: T, new_min: T, new_max: T) -> T
    where T: Copy + std::ops::Add<Output = T> + std::ops::Sub<Output = T> + std::ops::Div<Output = T>+ std::ops::Mul<Output = T>
{
    (((value - old_min) / (old_max - old_min)) * (new_max - new_min)) + new_min
}

// copied from fft_2d because it is not a public function:
/// Transposes a matrix represented as a flat 1D array in row-major order.
///
/// This function converts a matrix of size `height x width` (stored in a flat array) into its
/// transpose, returning the transposed matrix as a new `Vec<T>`.
///
/// # Parameters:
/// - `width`: Number of columns in the matrix.
/// - `height`: Number of rows in the matrix.
/// - `matrix`: The matrix data in row-major order.
///
/// # Returns:
/// A `Vec<T>` representing the transposed matrix.
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

/// Sums the diagonals of a 2D matrix (represented as a flat vector in row-major order).
///
/// This function computes the sum of the diagonals starting from the top-left corner and moving
/// to the bottom-right, for a given matrix with specified width and height.
///
/// # Parameters:
/// - `input`: A slice representing the matrix in row-major order.
/// - `width`: The number of columns in the matrix.
/// - `height`: The number of rows in the matrix.
///
/// # Returns:
/// A `Vec<T>` containing the sums of each diagonal.
///
/// # Example:
/// ```rust
/// let matrix = vec![1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12];
/// let result = sum_matrix_diagonals(&matrix, 4, 3);
/// assert_eq!(result, vec![1, 7, 18, 21, 19, 12]);
/// ```
pub fn sum_matrix_diagonals<T>(input: &[T], width: usize, height: usize) -> Vec<T>
    where T: Num + Copy + From<f64>
{
    let number_of_diagonals = width + (height - 1);
    let mut result: Vec<T> = Vec::with_capacity(number_of_diagonals);
    let mut offset: isize = -1;

    for (i, val) in input.iter().enumerate() {
        let i_row = i % width;
        if i_row == 0 { offset += 1 }
        let idx = i_row + offset as usize;

        let value = *val;

        if idx >= result.len() {
            result.push(value);
        } else {
            result[idx] = result[idx] + value;
        }
    }

    result
}

/// Expands a row-major linear slice into overlapping diagonal windows.
///
/// Each output row represents a diagonal slice of the input matrix, starting from each row index.
/// The diagonals maintain the original row order but shift right with increasing rows.
///
/// # Arguments
///
/// * `input` - A slice representing a matrix in row-major order.
/// * `width` - The number of columns in the matrix.
/// * `height` - The number of rows in the matrix.
///
/// # Returns
///
/// A `Vec<T>` containing overlapping diagonal windows extracted from the matrix.
///
/// # Example
///
/// ```
/// let input = vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0];
/// let result = linear_to_diagonals(&input, 4, 3);
/// assert_eq!(result, vec![1.0, 2.0, 3.0, 4.0,
///                         2.0, 3.0, 4.0, 5.0,
///                         3.0, 4.0, 5.0, 6.0]);
/// ```
pub fn linear_to_diagonals<T>(input: &[T], width: usize, height: usize) -> Vec<T>
    where T: Num + Copy + From<f64>
{
    if input.len() < ((width + height) - 1) {
        panic!("linear_to_diagonals: input too short for width and height!")
    }

    let mut result: Vec<T> = Vec::with_capacity(width * height);

    for (offset, _) in (0..height).enumerate() {
        for i in 0..width {
            result.push(input[i + offset])
        }
    }

    result

}

/// Performs downsampling in the frequency domain using FFT.
///
/// This function reduces the length of a complex signal by applying a forward FFT,
/// truncating high-frequency components, and then performing an inverse FFT.
/// It only downsamples when `new_length` is smaller than the original length,
/// preventing unintended upsampling.
///
/// # Arguments
///
/// * `input` - A slice of complex numbers representing the signal.
/// * `new_length` - The desired length after downsampling.
///
/// # Returns
///
/// A `Vec<Complex<f64>>` containing the downsampled signal.
///
/// # Notes
///
/// * The function applies a high-pass filter by truncating the frequency components.
/// * If `new_length >= input.len()`, the function returns a copy of the input.
///
/// # Example
///
/// ```
/// use rustfft::num_complex::Complex;
///
/// let signal = vec![Complex::new(1.0, 0.0), Complex::new(2.0, 0.0), Complex::new(3.0, 0.0), Complex::new(4.0, 0.0)];
/// let downsampled = downsampling(&signal, 2);
/// assert_eq!(downsampled.len(), 2);
/// ```
pub fn downsampling(input: &[Complex<f64>], new_length: usize) -> Vec<Complex<f64>> {
    // only do work when new_length < current_length (because we don't want to upsample)
    if new_length >= input.len() { return input.to_vec() }

    let current_length = input.len();
    let mut buffer: Vec<Complex<f64>> = input.to_vec();

    // FFT
    let mut planner_for = FftPlanner::new();
    let fft_for = planner_for.plan_fft_forward(current_length);
    fft_for.process(&mut buffer);

    // highpass
    buffer = buffer[0..new_length].to_vec();

    // inverse FFT
    let mut planner_inv = FftPlanner::new();
    let fft_inv = planner_inv.plan_fft_inverse(new_length);
    fft_inv.process(&mut buffer);

    buffer
}

// Same as above but for Vec<f64>
/*pub fn downsampling(input: &mut [f64], new_length: usize) -> Vec<f64> {
    // only do work when new_length < current_length (because we don't want to upsample)
    if new_length >= input.len() { return input.to_vec() }

    let current_length = input.len();

    // Convert the input buffer to complex numbers to be able to compute the FFT.
    let mut buffer: Vec<Complex<f64>> = input
        .iter()
        .map(|&val| Complex::new(val, 0.0))
        .collect();

    // FFT
    let mut planner_for = FftPlanner::new();
    let fft_for = planner_for.plan_fft_forward(current_length);
    fft_for.process(&mut *buffer);

    // highpass
    buffer = buffer[0..new_length].to_vec();

    // inverse FFT
    let mut planner_inv = FftPlanner::new();
    let fft_inv = planner_inv.plan_fft_inverse(new_length);
    fft_inv.process(&mut *buffer);

    get_magnitudes(&mut buffer)
}*/

/// Computes the magnitudes of complex numbers.
///
/// Calculates the Euclidean magnitude (absolute value) for each complex number in the input slice.
///
/// # Arguments
///
/// * `input` - A slice of complex numbers.
///
/// # Returns
///
/// A `Vec<f64>` containing the magnitudes of the input values.
///
/// # Example
///
/// ```
/// use rustfft::num_complex::Complex;
///
/// let input = vec![Complex::new(3.0, 4.0), Complex::new(1.0, 1.0)];
/// let magnitudes = get_magnitudes(&input);
/// assert_eq!(magnitudes, vec![5.0, (2.0f64).sqrt()]);
/// ```
pub fn get_magnitudes(input: &[Complex<f64>]) -> Vec<f64> {
    input.iter().map(|x| (pow(x.re, 2) + pow(x.im, 2)).sqrt()).collect()
}

/// Returns the maximum absolute value from a slice of `f64` values.
///
/// Computes the maximum of the absolute values in the input slice.
/// If the slice is empty, returns `0.0`.
pub fn max_abs(input: &[f64]) -> f64 {
    input.iter().fold(0.0_f64, |a, &b| a.max(b.abs()))
}

// This is currently not needed but let's keep for now...
/*
// take a 32x32 sub-image from an image buffer in row-major order.
fn partial_image_32<T: Copy>(image_data: &[T], width: usize, height: usize, (start_x, start_y): (usize, usize)) -> Result<Vec<T>, ()> {
    if start_x + 32 >= width || start_y + 32 >= height || width * height > image_data.len() {
       return Err(())
    }

    let mut result = Vec::with_capacity(32 * 32);

    // go through 32 rows
    for y in start_y..(start_y + 32) {
        let x = start_x + (y * width);
        for pixel in image_data.iter().skip(x).take(32) {
            result.push(*pixel)
        }
    };

    Ok(result)
}

// split one image vector into a collection of image vectors
pub fn split_image<T: Copy>(image_data: &[T], width: usize, height: usize) -> Vec<Vec<T>> {
    let mut result: Vec<Vec<T>> = Vec::new();
    let per_row = width / 32;
    let per_column = height / 32;

    for x in 0..per_row {
        for y in 0..per_column {
            result.push(partial_image_32(image_data, width, height, (x * 32, y * 32)).unwrap())
        }
    }

    result
}*/