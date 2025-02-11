use std::cmp::min;
use num_traits::Num;
use rustfft::FftPlanner;
use rustfft::num_complex::Complex;

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

pub fn _split_image<T: Copy>(image_data: &[T], width: usize, height: usize) -> Vec<Vec<T>> {
    let mut result: Vec<Vec<T>> = Vec::new();
    let per_row = width / 32;
    let per_column = height / 32;

    for x in 0..per_row {
        for y in 0..per_column {
            result.push(partial_image_32(image_data, width, height, (x * 32, y * 32)).unwrap())
        }
    }

    result
}

// Take a vector representing a 2D Matrix in row major order and sum all the diagonals, starting
// from the "top left" corner:
// sum_diagonals(vec![1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12], 4, 3) -> [1, 7, 18, 21, 19, 12]
pub fn sum_diagonals<T>(input: &[T], width: usize, height: usize) -> Vec<T>
    where T: Num + Copy + From<f64>
{
    let number_of_diagonals = width + (height - 1);
    // let nod_half: f64 = number_of_diagonals as f64 / 2.0;
    // let max_diagonal_size = min(width, height);
    let mut result: Vec<T> = Vec::with_capacity(number_of_diagonals);
    let mut offset: isize = -1;

    for (i, val) in input.iter().enumerate() {
        // position in result vector (index of diagonal)
        let i_row = i % width;
        if i_row == 0 { offset += 1 }
        let idx = i_row + offset as usize;

        // the following code would be used to average all values in a diagonal, but we don't need that
        /*
        // how many numbers are in this diagonal?
        // divide by this number to get the average after summing
        let nums_in_diagonal: T = T::from(
            if idx >= nod_half.floor() as usize {
                nod_half.round() as usize - ((idx + 1) % (nod_half.floor() as usize + 1))
            } else { 1 + idx }
            .min(max_diagonal_size) as f64);
            */

        let value = *val;  //    / nums_in_diagonal;

        if idx >= result.len() {
            result.push(value);
        } else {
            result[idx] = result[idx] + value;
        }
    }

    result
}

// If we also took the average, some tests for this would be:
// assert_eq!(vec![1.0, 3.5, 6.0, 7.0, 9.5, 12.0], average_diagonals(&vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0], 4, 3));
// assert_eq!(vec![1.0, 4.5, 5.5, 6.5, 7.5, 8.5, 12.0], average_diagonals(&vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0], 6, 2));
// assert_eq!(vec![1.0, 3.0, 5.0, 8.0, 10.0, 12.0], average_diagonals(&vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0], 3, 4));
// assert_eq!(vec![1.0, 4.0, 7.0, 8.0, 9.0, 12.0, 15.0], average_diagonals(&vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0], 5, 3));

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
    fft_inv.process(&mut *buffer);

    buffer
}