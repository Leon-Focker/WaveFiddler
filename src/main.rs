use std::cmp::min;
use rustfft::num_complex::Complex;
use rustfft::FftPlanner;
use fft2d::slice::fft_2d;
use plotters::prelude::*;
use rustfft::num_traits::pow;
use std::i16;
use hound;

fn main() -> Result<(), Box<dyn std::error::Error>>{

    // important parameters:
    let file_path = "data/hoch.png";
    let mut table_len = 0;
    let mut nr_harmonics = 0;

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
    // Now the matrix is organized as is visualized here: +
    // https://dsp.stackexchange.com/questions/3511/harmonics-in-2-d-fft-signal
    let transformed_img = transpose(width as usize, height as usize, &mut img_buffer);

    // Only take the harmonics we want
    // if table_len == 0 { table_len = transformed_img.len() }
    if table_len == 0 { table_len = transformed_img.len()}
    if nr_harmonics == 0 || nr_harmonics > table_len || nr_harmonics > transformed_img.len() {
        nr_harmonics = min(table_len, transformed_img.len())
    }
    // get the respective number of harmonics
    let mut buffer: Vec<Complex<f64>> = transformed_img[0..nr_harmonics].to_vec();

    // get magnitude of complex numbers and plot spectrogram
    let mag: Vec<f64> = buffer.clone().into_iter().map(|x| (pow(x.re, 2) + pow(x.im, 2)).sqrt()).collect();
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


fn _plot_complex_numbers(input: &Vec<Complex<f64>>) -> Result<(), Box<dyn std::error::Error>> {
    let x_min = input.iter().fold(f64::INFINITY, |a, &b| a.min(b.re));
    let x_max = input.iter().fold(f64::NEG_INFINITY, |a, &b| a.max(b.re));
    let y_min = input.iter().fold(f64::INFINITY, |a, &b| a.min(b.im));
    let y_max = input.iter().fold(f64::NEG_INFINITY, |a, &b| a.max(b.im));

    // Create a new chart
    let root = BitMapBackend::new("complex_plot.png", (640, 480)).into_drawing_area();
    root.fill(&WHITE)?;

    // Set up the chart
    let mut chart = ChartBuilder::on(&root)
        .caption("Simple Vector Plot", ("Arial", 40))
        .build_cartesian_2d(x_min..x_max, y_min.round()..y_max.round())?;

    // start with origin
    chart.draw_series(LineSeries::new(vec![(0.0, 0.0)].into_iter(), &BLUE).point_size(5))?;

    // Draw the x-y scatter plot
    chart.draw_series(LineSeries::new(
        input.into_iter().map(|complex| (complex.re, complex.im)),
        &RED,
    ))?;

    Ok(())
}

fn plot_numbers(name: &str, input: &Vec<f64>) -> Result<(), Box<dyn std::error::Error>> {
    let y_min = input.iter().fold(f64::INFINITY, |a, &b| a.min(b));
    let y_max = input.iter().fold(f64::NEG_INFINITY, |a, &b| a.max(b));

    // Create a new chart
    let root = BitMapBackend::new(name, (640, 480)).into_drawing_area();
    root.fill(&WHITE)?;

    // Set up the chart
    let mut chart = ChartBuilder::on(&root)
        .caption("Simple Vector Plot", ("Arial", 40))
        .build_cartesian_2d(0..(input.len() as i32), y_min.round()..y_max.round())?;

    // Draw the x-y scatter plot
    chart.draw_series(LineSeries::new(
        input.into_iter().enumerate().map(|(x, y)| (x as i32, *y)),
        &RED,
    ))?;

    Ok(())
}

fn write_to_wav (name: &str, data: &Vec<f64>) -> Result<(), Box<dyn std::error::Error>>  {

    let spec = hound::WavSpec {
        channels: 1,
        sample_rate: 48000,
        bits_per_sample: 16,
        sample_format: hound::SampleFormat::Int,
    };
    let mut writer = hound::WavWriter::create(name, spec).unwrap();

    let max = data.iter().fold(0.0, |a: f64, &b| a.max(b.abs()));
    let min = data.iter().fold(f64::INFINITY, |a: f64, &b| a.min(b.abs()));

    let rescaled = &data.iter().map(|x| rescale(*x, min, max, -(i16::MAX as f64), i16::MAX as f64)).collect();
    let _ = plot_numbers("wavetable.png", &rescaled);

    for sample in rescaled.into_iter() {
        writer.write_sample(*sample as i16).unwrap();
    }

    Ok(())
}

fn rescale<T: Copy + std::ops::Add<Output = T> + std::ops::Sub<Output = T> + std::ops::Div<Output = T>+ std::ops::Mul<Output = T>>
(value: T, old_min: T, old_max: T, new_min: T, new_max: T) -> T {
    (((value - old_min) / (old_max - old_min)) * (new_max - new_min)) + new_min
}

// copied from fft_2d because it is not a public function:
fn transpose<T: Copy + Default>(width: usize, height: usize, matrix: &[T]) -> Vec<T> {
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