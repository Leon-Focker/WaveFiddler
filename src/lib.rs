use rustfft::num_complex::Complex;
use rustfft::FftPlanner;
use fft2d::slice::{fft_2d, ifft_2d};
use crate::plot::plot_numbers;
use crate::from_to_wav::{get_wav_specs, read_from_wav, write_to_wav};
use crate::from_to_image::write_to_image;
use crate::utilities::*;
use std::path::Path;
use clap::Parser;

mod plot;
mod from_to_wav;
mod from_to_image;
mod utilities;

/// Command-line arguments for the application.
///
/// Defines options for input file, frame rate, FFT window size, and verbosity.
/// Additional options for image generation and audio processing are planned (TODO).
#[derive(Parser)]
#[command(version, about, long_about = None)]
pub struct Cli {
    /// Path to the input file
    pub input: Option<String>,

    /// The framerate that is used to generate an image sequence
    #[arg(short, long, default_value_t = 24)]
    frame_rate: u8,

    /// Window size of the fft
    #[arg(short = 's', long, default_value_t = 512)]
    fft_size: u32,

    // TODO which method for image generation
    // TODO select how many images or wavetables to generate
    // TODO decide how audio channels map to color
    // TODO normalisieren an/aus (für wavetables und bilder)

    /// Print information about what's happening
    #[arg(short, long, default_value_t = false)]
    verbose: bool,

    /// Plot the generated waveforms
    #[arg(long, default_value_t = false)]
    plot_waveforms: bool,

    /// Plot the spectra of the generated waveforms
    #[arg(long, default_value_t = false)]
    plot_spectra: bool,
}

//////////////////////////////////////////////////////////////////
//                    Sound to Image                            //
// ------------------------------------------------------------ //

// TODO (nicht) normalisieren...
/// Generate a sequence of images from the audio file data, specifics on how this data influences
/// the visuals can be supplied with the cli argument.
pub fn sound_to_img_sequence(file_path: &str, cl_arguments: &Cli) -> Result<(), Box<dyn std::error::Error>> {
    let file_name = Path::new(file_path).file_stem().and_then(|s| s.to_str()).unwrap();
    let samples = read_from_wav(file_path);
    let nr_samples = samples.len() as u32;
    let audio_spec = get_wav_specs(file_path);
    let sample_rate = audio_spec.sample_rate;
    let nr_channels = audio_spec.channels;
    let frame_rate = cl_arguments.frame_rate as u32;
    // we double table_size, so the windows overlap
    let table_size = (sample_rate / frame_rate) * 2;
    let nr_frames = ((nr_samples / nr_channels as u32) / (table_size / 2)) - 1;
    let width = (table_size - 1) / 2;
    let height= width;
    let verbose = cl_arguments.verbose;

    // Make sure a stereo sound file ist used!
    if nr_channels != 2 {
        return Err("please use an audio file with exactly 2 channels!".into());
    }

    // Generate all frames
    // TODO Select which ones to generate
    for i in 0..nr_frames {
        if verbose { println!("generating frame {i}"); }

        // dividing by 2 for overlapping windows and multiplying by 2 for 2 interleaved audio
        // channels cancels out here:
        let start_sample = (table_size * i) as usize;

        // ... here it does not, multiply by 2, for 2 audio channels
        let end_sample = start_sample + (table_size * 2) as usize;

        generate_frame_from_audio(
            format!("frames/{}_{}.jpg", file_name, i).as_str(),
            &samples[start_sample..end_sample],
            table_size as usize,
            width as usize,
            height as usize,
            cl_arguments)?;
    }

    Ok(())
}

/// Generate an image from the given audio data, which is assumed to be an interleaved stereo file
fn generate_frame_from_audio(
    name: &str,
    audio_data: &[f64],
    table_len: usize,
    width: usize,
    height: usize,
    _cl_arguments: &Cli)
    -> Result<(), Box<dyn std::error::Error>>
{
    // This will hold the colour values of all pixels
    // three elements in this Vec make one pixel in rgb (thus multiply by 3)
    let mut pixel_data: Vec<u8> = Vec::with_capacity(width * height * 3);

    // get all samples of the respective channels
    let left_channel: Vec<f64> = audio_data.iter().step_by(2).cloned().collect();
    let right_channel: Vec<f64> = audio_data.iter().skip(1).step_by(2).cloned().collect();
    let mid: Vec<f64> = audio_data.chunks(2)
        .map(|chunk| chunk.iter().sum())
        .collect();

    // TODO how do these generate a color, influenced by cl_arguments
    // TODO normalize or bold colours?
    let lum_vec1 = map_samples_into_2d(&left_channel, table_len);
    let lum_vec2 = map_samples_into_2d(&right_channel, table_len);
    let lum_vec3 = map_samples_into_2d(&mid, table_len);

    // Fill pixel_data with pixels
    // TODO see above
    for ((r, g), b) in lum_vec1.iter().zip(lum_vec2).zip(lum_vec3) {
        let r_val = (r * 255.0) as u8;
        let g_val = (g * 255.0) as u8;
        let b_val = (b * 255.0) as u8;
        //let w_val = (w / samples4_max * 255.0) as u8;
        pixel_data.push(r_val.saturating_add(0));
        pixel_data.push(g_val.saturating_add(0));
        pixel_data.push(b_val.saturating_add(0));
    }

    write_to_image(name, &pixel_data, width, height)?;

    Ok(())
}

/// Maps 1D samples into a 2D space using FFT and inverse FFT.
///
/// This function normalizes the input samples, performs a forward FFT (one-dimensional),
/// converts the result to a 2D form, applies an inverse FFT (two-dimensional), and returns the real
/// part of the resulting 2D samples.
///
/// # Arguments
///
/// * `samples` - A slice of 1D samples to be mapped.
/// * `table_len` - The desired length of the FFT table (used to determine 2D dimensions).
///
/// # Returns
///
/// A `Vec<f64>` representing the real part of the 2D FFT-transformed samples.
fn map_samples_into_2d(samples: &[f64], table_len: usize) -> Vec<f64> {
    let sample_max = max_abs(samples);
    let width = (table_len - 1) / 2;
    let height = width;

    // Convert the image buffer to complex numbers to be able to compute the FFT.
    let mut buffer: Vec<Complex<f64>> = samples
        .iter()
        .map(|&sample| Complex::new(sample / sample_max, 0.0))
        .collect();

    // apply the fft
    let mut planner = FftPlanner::new();
    let fft = planner.plan_fft_forward(table_len);
    fft.process(&mut buffer);

    // get 2D vector:
    // TODO other methods to get 2D from 1D?
    buffer = linear_to_diagonals(&buffer, width, height);

    // apply inverse 2D FFT
    ifft_2d(width, height, &mut buffer);

    // Output of ifft_2d is transposed, transpose back and collect real part
    transpose(width, height, &buffer)
        .iter()
        .map(|complex| complex.re )
        .collect()
}

//////////////////////////////////////////////////////////////////
//                    Image to Sound                            //
// ------------------------------------------------------------ //

/// Get rgb pixel data from an image and generate .wav files (wavetables)
pub fn image_to_waves(file_path: &str, cl_arguments: &Cli) -> Result<(), Box<dyn std::error::Error>> {
    let file_name = Path::new(file_path).file_stem().and_then(|s| s.to_str()).unwrap();
    let table_len = cl_arguments.fft_size.try_into().unwrap();

    // Open image from disk as rgb
    let img = image::open(file_path)?.into_rgb8();
    let (width, height) = img.dimensions();

    // Save the pixel_data into respective vectors
    let color_vec_len = img.len() / 3;
    let mut r_vector = Vec::with_capacity(color_vec_len);
    let mut g_vector = Vec::with_capacity(color_vec_len);
    let mut b_vector = Vec::with_capacity(color_vec_len);

    for pixel in img.as_raw().chunks_exact(3) {
        r_vector.push(pixel[0]);
        g_vector.push(pixel[1]);
        b_vector.push(pixel[2]);
    }

    // Generate the Wavetable
    // TODO decide via cl_arguments, whether to generate one wavetable frame or several
    // TODO Method to generate wavetables.
    let mut wavetable = Vec::with_capacity(table_len * 3);

    // get a wavetable (1D Vector of f64) from the image data (2D Vector of f64)
    wavetable.append(&mut waveform_from_image_data(r_vector, width as usize, height as usize, table_len, cl_arguments, 1));
    wavetable.append(&mut waveform_from_image_data(b_vector, width as usize, height as usize, table_len, cl_arguments, 2));
    wavetable.append(&mut waveform_from_image_data(g_vector, width as usize, height as usize, table_len, cl_arguments, 3));

    // write wavetable to audio file
    write_to_wav(
        format!("{}_{}.wav", file_name, table_len).as_str(),
        &wavetable,
        true,
        cl_arguments.plot_waveforms)?;

    Ok(())
}

/// Generates a waveform from image data using 2D FFT.
///
/// This function converts image pixel values to complex numbers, applies a 2D FFT,
/// normalizes the result, sums the diagonals, downsamples if necessary, and applies
/// an inverse FFT to generate a 1D waveform.
///
/// # Arguments
///
/// * `image_data` - A vector of pixel values (0-255) from the image.
/// * `width` - The width of the image.
/// * `height` - The height of the image.
/// * `table_len` - The desired length of the resulting waveform.
/// * `cl_arguments` - Command-line arguments for additional options (e.g., plotting).
/// * `id` - an ID to name the spectrum .png file.
///
/// # Returns
///
/// A vector of `f64` values representing the generated waveform.
fn waveform_from_image_data(
    image_data: Vec<u8>,
    width: usize,
    height: usize,
    mut table_len: usize,
    cl_arguments: &Cli,
    id: u8)
    -> Vec<f64> {

    if table_len == 0 { table_len = width + (height - 1)}
    let norm_factor: f64 = 1.0 / ((width as f64 * height as f64).sqrt()); // TODO necessary?

    // Convert the image buffer to complex numbers to be able to compute the FFT.
    let mut img_buffer: Vec<Complex<f64>> = image_data
        .iter()
        .map(|&pix| Complex::new(pix as f64 / 255.0, 0.0))
        .collect();

    // apply 2D FFT
    fft_2d(width, height, &mut img_buffer);

    // Normalize all values // TODO necessary?
    for num in img_buffer.iter_mut() {
        *num *= norm_factor;
    }

    // Output of fft_2d is transposed. But since we sum all the diagonals anyway,
    // we don't necessarily need transpose the matrix back. If we were to transpose it,
    // it would be arranged as visualized here:
    // https://dsp.stackexchange.com/questions/3511/harmonics-in-2-d-fft-signal
    /* let mut img_buffer = transpose(width, height, &img_buffer); */

    // Group all related harmonics together (sum_diagonals)
    // if necessary, reduce the buffer by downsampling.
    let mut buffer: Vec<Complex<f64>> = downsampling(
        &sum_matrix_diagonals(&img_buffer, width, height),
        table_len);

    // Maybe plot spectrum
    if cl_arguments.plot_spectra {
        // get magnitude of complex numbers and plot spectrogram
        let mag: Vec<f64> = get_magnitudes(&buffer);
        let name = &format!("spectrum_{}.png", id);
        if let Err(_err) = plot_numbers(name, &mag) {
            println!("couldn't plt spectrum: {name}");
        }
    }

    // If table_len is bigger than the provided buffer, pad the buffer with 0.0s
    if table_len > buffer.len() {
        for _ in 0..(table_len - buffer.len()) {
            buffer.push(Complex::new(0.0, 0.0));
        }
    }

    // apply the inverse fft (one-dimensional)
    let mut planner = FftPlanner::new();
    let fft = planner.plan_fft_inverse(table_len);
    fft.process(&mut buffer);

    // use only the real part as wavetable
    buffer.into_iter().map(|x| x.re).collect()
}