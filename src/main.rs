use std::cmp::min;
use rustfft::num_complex::Complex;
use rustfft::FftPlanner;
use fft2d::slice::fft_2d;
use rustfft::num_traits::pow;
use crate::plot::plot_numbers;
use crate::to_wav::write_to_wav;
use crate::utilities::*;

mod plot;
mod to_wav;
mod utilities;

fn main() -> Result<(), Box<dyn std::error::Error>>{

    // important parameters:
    let file_path = "data/mandrill.jpg";
    let table_len = 256;
    let nr_harmonics = 0;
    // 0 -> any number of harmonics, 1 -> the vertical only harmonics, 2 -> the horizontal only harmonics
    // see https://dsp.stackexchange.com/questions/3511/harmonics-in-2-d-fft-signal
    let method = 4;
    let normalize = true;

    // Open image from disk.
    let img = image::open(file_path)?.into_luma8();
    let (width, height) = img.dimensions();

    /*let test = split_image(img.clone().as_raw(), width as usize, height as usize);

    let mut real: Vec<f64> = Vec::new();

    for i in 0..test.len() {
        real.append(&mut transform_image_data(test[i].clone(), 32, 32, table_len, nr_harmonics, method));
    }*/

    let real = transform_image_data(img.as_raw().clone(), width, height, table_len, nr_harmonics, method);

    // write wavetable
    let _ = write_to_wav("wave.wav", &real, normalize);

    Ok(())
}

fn transform_image_data(image_data: Vec<u8>, width: u32, height: u32, table_len: usize, nr_harmonics: usize, method: i32) -> Vec<f64> {
    let norm_factor: f64 = 1.0 / ((width as f64 * height as f64).sqrt());
    // Convert the image buffer to complex numbers to be able to compute the FFT.
    let mut img_buffer: Vec<Complex<f64>> = image_data
        .iter()
        .map(|&pix| Complex::new(pix as f64 / 255.0, 0.0))
        .collect();

    // apply 2D FFT
    fft_2d(width as usize, height as usize, &mut img_buffer);

    // Output of fft_2d is transposed, transpose back.
    // Then the matrix is organized as is visualized here:
    // https://dsp.stackexchange.com/questions/3511/harmonics-in-2-d-fft-signal
    let mut transformed_img = transpose(width as usize, height as usize, &img_buffer);

    // Normalize all values
    for num in transformed_img.iter_mut() {
        *num *= norm_factor;
    }

    // Set DC offset to 0
    transformed_img[0] = Complex::new(0.0, 0.0);

    // get relevant harmonics and process them:
    match method {
        0 => generate_outputs_all(&transformed_img, nr_harmonics, table_len),
        1 => generate_outputs_vertical(&transformed_img, width as usize, table_len),
        2 => generate_outputs_horizontal(&transformed_img, width as usize, table_len),
        3 => chain_multiple_tables(&transformed_img, width as usize, height as usize, table_len),
        4 => generate_diagonal_outputs(&transformed_img, width as usize, height as usize, table_len),
        _ => panic!("invalid value for method")
    }
}

fn generate_diagonal_outputs(fft_output: &[Complex<f64>], width: usize, height: usize, mut table_len: usize) -> Vec<f64> {

    // adjust number of harmonics and table_len
    if table_len == 0 { table_len = width + (height - 1) }
    let nr_harmonics = min(table_len, width + (height - 1));

    //let mut buffer: Vec<Complex<f64>> = average_diagonals(fft_output, width, height)[0..nr_harmonics].to_vec();
    let mut buffer: Vec<Complex<f64>> = downsampling(&mut sum_diagonals(fft_output, width, height), table_len);

    // get magnitude of complex numbers and plot spectrogram
    let mag: Vec<f64> = get_magnitudes(&buffer);
    let _ = plot_numbers("spectrum.png", &mag);
    buffer = buffer[0..nr_harmonics].to_vec();

    // Pad the harmonic vector with Zeros, when we want fewer harmonics than the fft table is big
    if table_len > nr_harmonics {
        for _ in 0..(table_len - nr_harmonics) {
            buffer.push(Complex::new(0.0, 0.0));
        }
    }

    // apply the inverse fft
    let mut planner = FftPlanner::new();
    let fft = planner.plan_fft_inverse(table_len);
    fft.process(&mut buffer);

    // use only the real part as wave table
    buffer.into_iter().map(|x| x.re).collect()
}

fn chain_multiple_tables(fft_output: &[Complex<f64>], width: usize, height: usize, mut table_len: usize) -> Vec<f64> {

    // adjust number of harmonics and table_len
    if table_len == 0 { table_len = width }
    let nr_harmonics = min(table_len, width);

    let mut chained: Vec<f64> = Vec::new();

    for i in 0..height {
        let start_index = i * width;
        let end_index = (i + 1) * width;
        let mut buffer: Vec<Complex<f64>> = fft_output[start_index..end_index].to_vec();

        // Pad the harmonic vector with Zeros, when we want fewer harmonics than the fft table is big
        if table_len > nr_harmonics {
            for _ in 0..(table_len - nr_harmonics) {
                buffer.push(Complex::new(0.0, 0.0));
            }
        }

        // apply the inverse fft
        let mut planner = FftPlanner::new();
        let fft = planner.plan_fft_inverse(table_len);
        fft.process(&mut buffer);

        // use only the real part as wave table
        let mut real: Vec<f64> = buffer.into_iter().map(|x| x.re).collect();

        // rescale
        let max = real.iter().fold(0.0, |a: f64, &b| a.max(b.abs()));
        let min = real.iter().fold(f64::INFINITY, |a: f64, &b| a.min(b.abs()));

        // something is off with this rescaling...
        real = real.iter().map(|x| rescale(*x, min, max, 0.0, 1.0)).collect();

        chained.append(&mut real);
    }

    chained
}

fn generate_outputs_horizontal(fft_output: &[Complex<f64>], width: usize, mut table_len: usize) -> Vec<f64> {
    let mut buffer: Vec<Complex<f64>> = get_horizontal_harmonics(fft_output, width);

    // adjust number of harmonics and table_len
    if table_len == 0 { table_len = width }
    let nr_harmonics = min(table_len, buffer.len());

    // get magnitude of complex numbers and plot spectrogram
    let mag: Vec<f64> = get_magnitudes(&buffer);
    let _ = plot_numbers("spectrum.png", &mag);

    // Pad the harmonic vector with Zeros, when we want fewer harmonics than the fft table is big
    if table_len > nr_harmonics {
        for _ in 0..(table_len - nr_harmonics) {
            buffer.push(Complex::new(0.0, 0.0));
        }
    }

    // apply the inverse fft
    let mut planner = FftPlanner::new();
    let fft = planner.plan_fft_inverse(table_len);
    fft.process(&mut buffer);

    // use only the real part as wave table
    buffer.into_iter().map(|x| x.re).collect()
}

fn generate_outputs_vertical(fft_output: &[Complex<f64>], width: usize, mut table_len: usize) -> Vec<f64>  {
    let mut buffer: Vec<Complex<f64>> = get_vertical_harmonics(fft_output, width);

    // adjust number of harmonics and table_len
    if table_len == 0 { table_len = width }
    let nr_harmonics = min(table_len, buffer.len());

    // get magnitude of complex numbers and plot spectrogram
    let mag: Vec<f64> = get_magnitudes(&buffer);
    let _ = plot_numbers("spectrum.png", &mag);

    // Pad the harmonic vector with Zeros, when we want fewer harmonics than the fft table is big
    if table_len > nr_harmonics {
        for _ in 0..(table_len - nr_harmonics) {
            buffer.push(Complex::new(0.0, 0.0));
        }
    }

    // apply the inverse fft
    let mut planner = FftPlanner::new();
    let fft = planner.plan_fft_inverse(table_len);
    fft.process(&mut buffer);

    // use only the real part as wave table
    buffer.into_iter().map(|x| x.re).collect()
}

fn generate_outputs_all(fft_output: &[Complex<f64>], mut nr_harmonics: usize, mut table_len: usize) -> Vec<f64>  {

    // adjust number of harmonics and table_len
    if table_len == 0 { table_len = fft_output.len() }
    if nr_harmonics == 0 || nr_harmonics > table_len || nr_harmonics > fft_output.len() {
        nr_harmonics = min(table_len, fft_output.len())
    }
        
    let mut buffer: Vec<Complex<f64>> = get_first_harmonics(fft_output, nr_harmonics);

    // get magnitude of complex numbers and plot spectrogram
    let mag: Vec<f64> = get_magnitudes(&buffer);
    let _ = plot_numbers("spectrum.png", &mag);

    // Pad the harmonic vector with Zeros, when we want fewer harmonics than the fft table is big
    if table_len > nr_harmonics {
        for _ in 0..(table_len - nr_harmonics) {
            buffer.push(Complex::new(0.0, 0.0));
        }
    }

    // apply the inverse fft
    let mut planner = FftPlanner::new();
    let fft = planner.plan_fft_inverse(table_len);
    fft.process(&mut buffer);

    // use only the real part as wave table
    buffer.into_iter().map(|x| x.re).collect()
}

fn get_magnitudes(input: &[Complex<f64>]) -> Vec<f64> {
    input.iter().map(|x| (pow(x.re, 2) + pow(x.im, 2)).sqrt()).collect()
}

fn get_horizontal_harmonics(input: &[Complex<f64>], width: usize) -> Vec<Complex<f64>> {
    input.iter()
        .enumerate()
        .filter(|(i, _)| i % width == 0)
        .map(|(_, elem)| *elem)
        .collect()
}

fn get_vertical_harmonics(input: &[Complex<f64>], width: usize) -> Vec<Complex<f64>> {
    input[0..width].to_vec()
}

fn get_first_harmonics(input: &[Complex<f64>], how_many: usize) -> Vec<Complex<f64>> {
    input[0..how_many].to_vec()
}