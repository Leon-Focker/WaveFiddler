use rustfft::num_complex::Complex;
use rustfft::FftPlanner;
use fft2d::slice::fft_2d;
use crate::plot::plot_numbers;
use crate::to_wav::write_to_wav;
use crate::utilities::*;
use std::path::Path;

mod plot;
mod to_wav;
mod utilities;

fn main() -> Result<(), Box<dyn std::error::Error>>{

    // important parameters:
    let file_path = "data/banana.jpg";
    let table_len = 512;

    // Open image from disk as rgb, save the color values into vectors
    let img = image::open(file_path)?.into_rgb8();
    let (width, height) = img.dimensions();
    let color_vec_len = img.len() / 3;
    let mut r_vector = Vec::with_capacity(color_vec_len);
    let mut g_vector = Vec::with_capacity(color_vec_len);
    let mut b_vector = Vec::with_capacity(color_vec_len);

    for pixel in img.as_raw().chunks_exact(3) {
        r_vector.push(pixel[0]);
        g_vector.push(pixel[1]);
        b_vector.push(pixel[2]);
    }

    let mut wavetable = Vec::with_capacity(table_len * 3);

    // get a wavetable (1D Vector of f64) from the image data (2D Vector of f64)
    wavetable.append(&mut get_wave_table_from_image(r_vector, width, height, table_len));
    wavetable.append(&mut get_wave_table_from_image(b_vector, width, height, table_len));
    wavetable.append(&mut get_wave_table_from_image(g_vector, width, height, table_len));

    // write wavetable to audio file
    let file_name = Path::new(file_path).file_stem().unwrap().to_str().unwrap();
    let _ = write_to_wav(format!("{}_{}{}", file_name, table_len, ".wav").as_str(), &wavetable, true);

    Ok(())
}

fn get_wave_table_from_image(image_data: Vec<u8>, width: u32, height: u32, mut table_len: usize) -> Vec<f64> {

    if table_len == 0 { table_len = (width + (height - 1)) as usize }
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

    // get relevant harmonics and process them:
    generate_diagonal_outputs(&transformed_img, width as usize, height as usize, table_len)
}

fn generate_diagonal_outputs(fft_output: &[Complex<f64>], width: usize, height: usize, table_len: usize) -> Vec<f64> {

    // Group all related harmonics together (sum_diagonals) and, if necessary, reduce the
    // buffer by downsampling.
    let mut buffer: Vec<Complex<f64>> =
        downsampling(&sum_diagonals(fft_output, width, height), table_len);

    // get magnitude of complex numbers and plot spectrogram
    let mag: Vec<f64> = get_magnitudes(&buffer);
    let _ = plot_numbers("spectrum.png", &mag);

    // If table_len is bigger than the provided buffer, pad the buffer with 0.0s
    if table_len > buffer.len() {
        for _ in 0..(table_len - buffer.len()) {
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


/*fn chain_multiple_tables(fft_output: &[Complex<f64>], width: usize, height: usize, mut table_len: usize) -> Vec<f64> {

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
*/