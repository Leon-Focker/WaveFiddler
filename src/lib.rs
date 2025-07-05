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
pub mod utilities;

/// Command-line arguments for the application.
#[derive(Parser)]
#[command(version, about, long_about = None)]
pub struct Cli {
    /// Path to the input file
    pub input: Option<String>,

    /// Path to the output directory. If it does not exist, it will be created.
    #[arg(short, long, default_value_t = String::from("./wavefiddler_outputs/"))]
    pub output_dir: String,

    /// Name of the output file. Valid file extensions are .wav, .png and .jpg
    #[arg(short, long)]
    pub name: Option<String>,

    /// Window size of the fft. When generating a wavetable, this will be the length of the wavetable. If it is 0,
    /// the fft-size is automatically set so that no downsampling or padding has to be applied (not necessarily a power of 2).
    /// When converting audio to an image, the fft-size is only relevant when generating a single frame,
    /// else it is set automatically. Setting it to 0 uses the entire .wav file (could be too long!).
    #[arg(short = 'f', long, default_value_t = 512)]
    fft_size: u32,

    /// For image->audio conversion. When true, try to retain the spectral envelope. This might introduce gaps into the resulting audio file.
    /// Else use the spectrum as is, and discard harmonics, if there is too many.
    #[arg(short = 'S', long, default_value_t = false)]
    stretch_spectrum: bool,

    /// Image->audio conversion method. 0 => generate one frame. 1 => generate three frames from the rgb channels.
    /// 2 or higher => generate this many frames from a subset of the image spectrum.
    #[arg(short = 'i', long, default_value_t = 0)]
    pub i2a_method: u64,

    /// The framerate that is used to generate an image sequence
    #[arg(short = 'F', long, default_value_t = 24)]
    frame_rate: u8,

    /// Generate a single frame instead of the entire image sequence.
    /// Specify the start sample with this argument. As many samples as defined by --fft_size wil be used.
    #[arg(short = 's', long)]
    pub single_frame: Option<u64>,

    /// Audio->image conversion method. 0 => simple diagonal map, 1 => diagonal map + inverse 2D fft,
    /// 2 => simple linear map (significantly smaller images with default dimensions!), 3 => linear map + inverse 2D fft,
    /// 4 => same as 2, but transposed (flip the matrix but don't swap width / height), 5 or higher => same as 3 but transposed.
    #[arg(short = 'a', long, default_value_t = 0)]
    pub a2i_method: u64,

    /// How the left and right channel of the .wav file influence the color of the generated image.
    /// This argument must be a list of 6 numbers, which create color like this:
    /// red = left * nr1 + right * nr2, green = left * nr3 + right * nr4, blue = left * nr5 + right * nr6
    #[arg(short, long, default_value_t = String::from("1,0,0,1,0.5,0.5"))]
    color_map: String,

    /// The dimensions for the generated image as a list of width,height.
    /// If none are specified, they are chosen automatically to produce a square image.
    /// This does not simple 'rescale' the generated image! Changing the image dimensions has
    /// a huge impact in the generation of the images. While method 2 might be less interesting,
    /// methods 2-5 yield way better results.
    #[arg(short = 'd', long)]
    image_dimensions: Option<String>,

    /// How vibrant a generated image will be. Values between 0 and 5 work best.
    #[arg(short = 'v', long, default_value_t = 3.0)]
    vibrancy: f64,

    /// Print information about what's happening
    #[arg(long, default_value_t = false)]
    verbose: bool,

    /// For image->audio conversion. Plot the generated waveforms
    #[arg(long, default_value_t = false)]
    plot_waveforms: bool,

    /// For image->audio conversion. Plot the spectra of the generated waveforms
    #[arg(long, default_value_t = false)]
    plot_spectra: bool,

}

//////////////////////////////////////////////////////////////////
//                    Sound to Image                            //
// ------------------------------------------------------------ //

/// Generate a sequence of images from the audio file data, specifics on how this data influences
/// the visuals can be supplied with the cli argument.
pub fn sound_to_img_sequence(file_path: &str, cl_arguments: &Cli) -> Result<(), Box<dyn std::error::Error>> {
    let file_name = Path::new(file_path).file_stem().and_then(|s| s.to_str()).unwrap();
    let samples = read_from_wav(file_path)?;
    let nr_samples = samples.len() as u32;
    let audio_spec = get_wav_specs(file_path);
    let sample_rate = audio_spec.sample_rate;
    let nr_channels = audio_spec.channels;
    let frame_rate = cl_arguments.frame_rate as u32;

    // we double table_size, so the windows overlap
    let table_size = (sample_rate / frame_rate) * 2;
    let nr_frames = ((nr_samples / nr_channels as u32) / (table_size / 2)).saturating_sub(1);

    // make sure the output directory exists
    ensure_directory_exists(&cl_arguments.output_dir)
        .expect("something went wrong while checking for the output directory");

    // Get output file name
    let out_file_name: &str = match &cl_arguments.name {
        None => &format!("{}.jpg", file_name),
        Some(string) => {
            let ext = Path::new(string).extension().and_then(|s| s.to_str()).unwrap_or("jpg");
            if !["jpg", "png"].contains(&ext) {
            return Err("invalid output file extension when converting to an image!".into());
            } else { &format!("{}.{}", Path::new(string).file_stem().and_then(|s| s.to_str()).unwrap(), ext) }
        },
    };

    // Make sure a stereo sound file ist used!
    if nr_channels > 2 {
        return Err("please use an audio file with not more than 2 channels!".into());
    }

    // Generate all frames
    match cl_arguments.single_frame {
        Some(sample_idx) => {
            // Information for the user
            println!("Generating an image from audio...");

            // multiply by 2 (stereo file)
            let start_sample = (sample_idx * 2) as usize;
            let mut table_size = cl_arguments.fft_size as usize;
            if table_size == 0 { table_size = nr_samples as usize / 2 };
            let end_sample = start_sample + (table_size * 2) - 1;

            if end_sample >= nr_samples as usize {
                return Err(format!("Starting at sample {} with an fft-size of {} is out of bounds for a file with {} samples!",
                start_sample, table_size, nr_samples)
                    .into());
            }

            generate_frame_from_audio(
                format!("{}{}", &cl_arguments.output_dir, out_file_name).as_str(),
                &samples[start_sample..end_sample],
                nr_channels,
                cl_arguments)?;
        },
        None => {
            // Information for the User
            println!("Generating an image sequence from audio...");
            println!(); // This will immediately be replaced by the progress bar

            for i in 0..nr_frames {
                // Progress information
                if print_progress_bar(
                    i as f32,
                    (nr_frames - 1) as f32,
                    50.0,
                    cl_arguments.verbose).is_err() {
                    println!("Printing progress bar failed!");
                }

                // dividing by 2 for overlapping windows and multiplying by nr_channels (interleaved audio)
                let start_sample = (table_size * i / 2) as usize * nr_channels as usize;
                let end_sample = start_sample + (table_size as usize * nr_channels as usize);

                generate_frame_from_audio(
                    format!("{}{:05}_{}", &cl_arguments.output_dir, i, out_file_name).as_str(),
                    &samples[start_sample..end_sample],
                    nr_channels,
                    cl_arguments)?;
            }
        },
    }

    Ok(())
}

/// Generate an image from the given audio data, which is assumed to be an interleaved stereo file
fn generate_frame_from_audio(
    name: &str,
    audio_data: &[f64],
    audio_channels: u16,
    cl_arguments: &Cli)
    -> Result<(), Box<dyn std::error::Error>>
{
    // get all samples of the respective channels
    let (left_channel, right_channel): (Vec<f64>, Vec<f64>) =
        if audio_channels == 1 {
            (audio_data.to_vec(), audio_data.to_vec())
        } else {
            (
                audio_data.iter().step_by(2).cloned().collect(),
                audio_data.iter().skip(1).step_by(2).cloned().collect(),
            )
        };

    // get the image dimensions
    let (width, height, resample): (usize, usize, bool) =
        match &cl_arguments.image_dimensions {
            // calculate width and height
            None => {
                if cl_arguments.a2i_method < 2 {
                    ((left_channel.len() - 1) / 2, (left_channel.len() - 1) / 2, false)
                } else {
                    (left_channel.len().isqrt(), left_channel.len() / left_channel.len().isqrt(), false)
                }
            }
            // use width, height from cli
            Some(s)=> {
                // from user argument
                let mut dimensions = s
                    .split(',')
                    .map(|x| x.parse::<usize>().unwrap_or(1));
                (dimensions.next().unwrap_or(1), dimensions.next().unwrap_or(1), true)
            }
        };

    // map those samples into 2D space
    let lum_vec_left = map_samples_into_2d(&left_channel, width, height, resample, cl_arguments);
    let lum_vec_right = map_samples_into_2d(&right_channel, width, height, resample, cl_arguments);
    // max value
    let lum_max = max_abs(&lum_vec_left).max(max_abs(&lum_vec_right))
        / f64::powf(10.0, cl_arguments.vibrancy);
    let scale = 255.0 / lum_max;

    // This will hold the colour values of all pixels
    // three elements in this Vec make one pixel in rgb (thus multiply by 3)
    let mut pixel_data: Vec<u8> = Vec::with_capacity(width * height * 3);

    // How the audio channels influence the color
    let color_mult: [f64; 6] = cl_arguments
        .color_map
        .split(',')
        .map(|x| x.parse::<f64>().unwrap_or(0.0))
        .collect::<Vec<f64>>()
        .try_into()
        .expect("Expected exactly 6 color multipliers");

    // Fill pixel_data with pixels
    for (&left, &right) in lum_vec_left.iter().zip(&lum_vec_right) {
        let l = left * scale;
        let r = right * scale;
        pixel_data.push(((l * color_mult[0]) + (r * color_mult[1])) as u8);
        pixel_data.push(((l * color_mult[2]) + (r * color_mult[3])) as u8);
        pixel_data.push(((l * color_mult[4]) + (r * color_mult[5])) as u8);
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
/// * `cl_arguments` - Command-line arguments for additional options (e.g., method selection).
///
/// # Returns
///
/// A `Vec<f64>` representing the real part of the 2D FFT-transformed samples.
fn map_samples_into_2d(samples: &[f64], width: usize, height: usize, needs_resample: bool, cl_arguments: &Cli) -> Vec<f64> {
    let sample_max = max_abs(samples);
    let table_len = samples.len();

    // Convert the image buffer to complex numbers to be able to compute the FFT.
    let mut buffer: Vec<Complex<f64>> = samples
        .iter()
        .map(|&sample| Complex::new(sample / sample_max, 0.0))
        .collect();

    // apply the fft
    let mut planner = FftPlanner::new();
    let fft = planner.plan_fft_forward(table_len);
    fft.process(&mut buffer);

    // resample when necessary
    if needs_resample {
        buffer = resample(&buffer, width * height);
    }

    if cl_arguments.a2i_method < 2 {
        // get 2D vector:
        buffer = linear_to_diagonals(&buffer, width, height);

        if cl_arguments.a2i_method == 1 {
            // apply inverse 2D FFT
            ifft_2d(width, height, &mut buffer);
        }
        buffer.iter().map(|complex| complex.re).collect()

    } else {

        buffer.truncate(width * height);

        if cl_arguments.a2i_method == 3 || cl_arguments.a2i_method == 5 {
            // apply inverse 2D FFT
            ifft_2d(width, height, &mut buffer);
        }

        if cl_arguments.a2i_method > 3 {
            // transpose the matrix
            transpose(width, height, &buffer)
                .iter()
                .map(|complex| complex.re)
                .collect()
        } else {
            buffer.iter().map(|complex| complex.re).collect()
        }
    }
}

//////////////////////////////////////////////////////////////////
//                    Image to Sound                            //
// ------------------------------------------------------------ //

/// Get rgb pixel data from an image and generate .wav files (wavetables)
pub fn image_to_waves(file_path: &str, cl_arguments: &Cli) -> Result<(), Box<dyn std::error::Error>> {
    let file_name = Path::new(file_path).file_stem().and_then(|s| s.to_str()).unwrap();
    let table_len = cl_arguments.fft_size.try_into().unwrap();
    let audio_data_len = match cl_arguments.i2a_method {
        0 => table_len,
        1 => table_len * 3,
        _ => table_len * cl_arguments.i2a_method as usize,
    };

    // make sure the output directory exists
    ensure_directory_exists(&cl_arguments.output_dir)
        .expect("something went wrong while checking for the output directory");

    // store the generated audio data here
    let mut audio_data = Vec::with_capacity(audio_data_len);

    // Get output file name
    let out_file_name: &str = match &cl_arguments.name {
        None => &format!("{}_{}.wav", file_name, table_len),
        Some(string) => {
            let ext = Path::new(string).extension().and_then(|s| s.to_str()).unwrap_or("wav");
            if !["wav"].contains(&ext) {
                return Err("invalid output file extension when converting to a audio!".into());
            } else { &format!("{}.{}", Path::new(string).file_stem().and_then(|s| s.to_str()).unwrap(), ext) }
        },
    };

    // Different methods for generating the audio data
    if cl_arguments.i2a_method == 1 {
        // Information for the User
        println!("Generating audio from the rgb data of an image...");

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

        // get a wavetable (1D Vector of f64) from the image data (2D Vector of f64) for each rgb channel
        audio_data.append(&mut waveform_from_image_data(r_vector, width as usize, height as usize, table_len, cl_arguments, "r"));
        audio_data.append(&mut waveform_from_image_data(b_vector, width as usize, height as usize, table_len, cl_arguments, "g"));
        audio_data.append(&mut waveform_from_image_data(g_vector, width as usize, height as usize, table_len, cl_arguments, "b"));

    } else {
        // Information for the User
        println!("audio from a greyscale image...");

        // Open image from disk as greyscale
        let img = image::open(file_path)?.into_luma8();
        let (width, height) = img.dimensions();
        let pixel_data = img.as_raw().clone();

        if cl_arguments.i2a_method == 0 {
            // get a wavetable (1D Vector of f64) from the image data (2D Vector of f64) for grey scale image
            audio_data = waveform_from_image_data(pixel_data, width as usize, height as usize, table_len, cl_arguments, "grey")
        } else {
            // get a multiple wavetables (1D Vector of f64) from the image data (2D Vector of f64) for grey scale image
            audio_data = multiple_waveforms_from_image_data(pixel_data, width as usize, height as usize, table_len, cl_arguments, cl_arguments.i2a_method as usize)
        }
    }

    // write audio_data to audio file
    write_to_wav(
        format!("{}{}", &cl_arguments.output_dir, out_file_name).as_str(),
        &audio_data,
        true,
        if cl_arguments.plot_waveforms {
            Some(format!("{}waveform.png", cl_arguments.output_dir))
        } else { None }
    )?;

    Ok(())
}

/// Generates a waveform from image data using 2D FFT.
///
/// This function converts image pixel values to complex numbers, applies a 2D FFT,
/// sums the diagonals, resamples if necessary, and applies
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
    id: &str)
    -> Vec<f64> {

    if table_len == 0 { table_len = width + (height - 1)}

    // Convert the image buffer to complex numbers to be able to compute the FFT.
    let mut img_buffer: Vec<Complex<f64>> = image_data
        .iter()
        .map(|&pix| Complex::new(pix as f64 / 255.0, 0.0))
        .collect();

    // apply 2D FFT
    fft_2d(width, height, &mut img_buffer);

    // remove DC offset
    img_buffer[0] = Complex::new(0.0, 0.0);

    // Output of fft_2d is transposed. But since we sum all the diagonals anyway,
    // we don't necessarily need transpose the matrix back. If we were to transpose it,
    // it would be arranged as visualized here:
    // https://dsp.stackexchange.com/questions/3511/harmonics-in-2-d-fft-signal
    /* let mut img_buffer = transpose(width, height, &img_buffer); */

    let mut buffer: Vec<Complex<f64>> = sum_matrix_diagonals(&img_buffer, width, height);

    // Different approaches to handling different buffer length / table length
    if cl_arguments.stretch_spectrum {
        // Stretch the buffer, so it fits the table length
        buffer = resample(&buffer, table_len);
    } else if buffer.len() > table_len {
        // If the buffer holds more elements than the table needs, discard them
        buffer.truncate(table_len);
    } else {
        // If table_len is bigger than the provided buffer, pad the buffer with 0.0s
        for _ in 0..(table_len - buffer.len()) {
            buffer.push(Complex::new(0.0, 0.0));

        }
    }

    // Maybe plot spectrum
    if cl_arguments.plot_spectra {
        // get magnitude of complex numbers and plot spectrogram
        let mag: Vec<f64> = get_magnitudes(&buffer);
        let name = &format!("{}spectrum_{}.png", &cl_arguments.output_dir, id);
        if let Err(_err) = plot_numbers(name, &mag) {
            println!("couldn't plot spectrum: {name}");
        }
    }

    // apply the inverse fft (one-dimensional)
    let mut planner = FftPlanner::new();
    let fft = planner.plan_fft_inverse(table_len);
    fft.process(&mut buffer);

    // use only the real part as wavetable
    buffer.into_iter().map(|x| x.re).collect()
}

/// Generates multiple waveforms from image data using 2D FFT.
///
/// This function converts image pixel values to complex numbers, applies a 2D FFT,
/// removes the DC offset, and splits the result into a number of diagonals through the 2D matrix.
/// Apply a 1D inverse FFT to each of those diagonals and append them to one audio_data Vector.
///
/// # Arguments
///
/// * `image_data` - A vector of pixel values (0-255) from the image.
/// * `width` - The width of the image.
/// * `height` - The height of the image.
/// * `table_len` - The desired length of the resulting waveform.
/// * `cl_arguments` - Command-line arguments for additional options (e.g., plotting).
/// * `number_of_frames` - The number of waveforms to generate.
///
/// # Returns
///
/// A vector of `f64` values representing the generated waveform. Length = table_size * number_of_frames.
fn multiple_waveforms_from_image_data(
    image_data: Vec<u8>,
    width: usize,
    height: usize,
    mut table_len: usize,
    cl_arguments: &Cli,
    number_of_frames: usize)
    -> Vec<f64> {

    if table_len == 0 { table_len = width + (height - 1)}

    // Convert the image buffer to complex numbers to be able to compute the FFT.
    let mut img_buffer: Vec<Complex<f64>> = image_data
        .iter()
        .map(|&pix| Complex::new(pix as f64 / 255.0, 0.0))
        .collect();

    // apply 2D FFT
    fft_2d(width, height, &mut img_buffer);

    // remove DC offset
    img_buffer[0] = Complex::new(0.0, 0.0);

    // Output of fft_2d is transposed.Transpose it back to be arranged as depicted here:
    // https://dsp.stackexchange.com/questions/3511/harmonics-in-2-d-fft-signal
    img_buffer = transpose(width, height, &img_buffer);

    let mut result: Vec<f64> = Vec::with_capacity(table_len * number_of_frames);

    // panic, because this call should never produce an error caused by user input...
    let mut buffer_vec: Vec<Vec<Complex<f64>>> =
        get_matrix_rays(&img_buffer, width, height, number_of_frames)
        .unwrap_or_else(|_err| {
            panic!("getting matrix rays went wrong!");
        });

    // plan the inverse fft
    let mut planner = FftPlanner::new();
    let fft = planner.plan_fft_inverse(table_len);

    println!(); // This will immediately be replaced by the progress bar

    // apply an inverse fft to all buffers in buffer_vec
    for (i, buffer) in buffer_vec.iter_mut().enumerate() {
        // Progress information
        if print_progress_bar(
            i as f32,
            (number_of_frames - 1) as f32,
            50.0,
            cl_arguments.verbose).is_err() {
            println!("Printing progress bar failed!");
        }

        // Different approaches to handling different buffer length / table length
        if cl_arguments.stretch_spectrum {
            // Stretch the buffer, so it fits the table length
            *buffer = resample(buffer, table_len);
        } else if buffer.len() > table_len {
            // If the buffer holds more elements than the table needs, discard them
            buffer.truncate(table_len);
        } else {
            // If table_len is bigger than the provided buffer, pad the buffer with 0.0s
            for _ in 0..(table_len - buffer.len()) {
                buffer.push(Complex::new(0.0, 0.0));
            }
        }

        // apply the inverse fft (one-dimensional)
        fft.process(buffer);

        // use only the real part as wavetable
        result.append(&mut buffer.iter().map(|x| x.re).collect::<Vec<f64>>());

    }

    result
}