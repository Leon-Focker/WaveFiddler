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
    let file_path = "data/hoch.png";
    let table_len = 0;
    let nr_harmonics = 0;
    // 0 -> any number of harmonics, 1 -> the vertical only harmonics, 2 -> the horizontal only harmonics
    // see https://dsp.stackexchange.com/questions/3511/harmonics-in-2-d-fft-signal
    let method = 1;

    // Open image from disk.
    let img = image::open(file_path)?.into_luma8();
    let (width, height) = img.dimensions();

    // Convert the image buffer to complex numbers to be able to compute the FFT.
    let mut img_buffer: Vec<Complex<f64>> = img
        .as_raw()
        .iter()
        .map(|&pix| Complex::new(pix as f64 / 255.0, 0.0))
        .collect();

    // apply 2D FFT
    fft_2d(width as usize, height as usize, &mut img_buffer);

    // Output of fft_2d is transposed, transpose back.
    // Now the matrix is organized as is visualized here:
    // https://dsp.stackexchange.com/questions/3511/harmonics-in-2-d-fft-signal
    let transformed_img = transpose(width as usize, height as usize, &mut img_buffer);
    
    // get relevant harmonics and process them:
    match method {
        0 => generate_outputs_all(&transformed_img, nr_harmonics, table_len),
        1 => generate_outputs_vertical(&transformed_img, width as usize, table_len),
        2 => generate_outputs_horizontal(&transformed_img, width as usize, table_len),
        _ => Err(Box::from("invalid value for method"))
    }
}


fn generate_outputs_horizontal(fft_output: &Vec<Complex<f64>>, width: usize, mut table_len: usize)
    -> Result<(), Box<dyn std::error::Error>> {
    let mut buffer: Vec<Complex<f64>> = get_horizontal_harmonics(&fft_output, width);

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
    let real: Vec<f64> = buffer.into_iter().map(|x| x.re).collect();

    let _ = write_to_wav("wave.wav", &real);

    Ok(())
}

fn generate_outputs_vertical(fft_output: &Vec<Complex<f64>>, width: usize, mut table_len: usize)
    -> Result<(), Box<dyn std::error::Error>> {
    let mut buffer: Vec<Complex<f64>> = get_vertical_harmonics(&fft_output, width);

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
    let real: Vec<f64> = buffer.into_iter().map(|x| x.re).collect();

    let _ = write_to_wav("wave.wav", &real);

    Ok(())
}

fn generate_outputs_all(fft_output: &Vec<Complex<f64>>, mut nr_harmonics: usize, mut table_len: usize)
    -> Result<(), Box<dyn std::error::Error>> {

    // adjust number of harmonics and table_len
    if table_len == 0 { table_len = fft_output.len() }
    if nr_harmonics == 0 || nr_harmonics > table_len || nr_harmonics > fft_output.len() {
        nr_harmonics = min(table_len, fft_output.len())
    }
        
    let mut buffer: Vec<Complex<f64>> = get_first_harmonics(&fft_output, nr_harmonics);

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
    let real: Vec<f64> = buffer.into_iter().map(|x| x.re).collect();

    let _ = write_to_wav("wave.wav", &real);

    Ok(())
}

fn get_magnitudes(input: &Vec<Complex<f64>>) -> Vec<f64> {
    input.clone().into_iter().map(|x| (pow(x.re, 2) + pow(x.im, 2)).sqrt()).collect()
}

fn get_horizontal_harmonics(input: &Vec<Complex<f64>>, width: usize) -> Vec<Complex<f64>> {
    input.iter()
        .enumerate()
        .filter(|(i, _)| i % width == 0)
        .map(|(_, elem)| *elem)
        .collect()
}

fn get_vertical_harmonics(input: &Vec<Complex<f64>>, width: usize) -> Vec<Complex<f64>> {
    input[0..width].to_vec()
}

fn get_first_harmonics(input: &Vec<Complex<f64>>, how_many: usize) -> Vec<Complex<f64>> {
    input[0..how_many].to_vec()
}