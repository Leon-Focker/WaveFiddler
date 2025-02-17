use hound::WavSpec;
use crate::utilities::{max_abs};
use crate::plot::plot_numbers;

/// Write numbers from a vector into a .wav file. Numbers should be between -1.0 and 1.0,
/// else, when normalize is true, they are rescaled to fit that range.
pub fn write_to_wav (name: &str, data: &[f64], normalize: bool, plot_waveform: Option<String>) -> Result<(), Box<dyn std::error::Error>>  {

    let spec = WavSpec {
        channels: 1,
        sample_rate: 48000,
        bits_per_sample: 16,
        sample_format: hound::SampleFormat::Int,
    };

    let mut writer = hound::WavWriter::create(name, spec).unwrap();

    // normalize if necessary
    let samples: Vec<f64> =  if normalize {
        let y_max = max_abs(&data);

        data.iter().map(|x| (*x / y_max) * i16::MAX as f64).collect()
    } else {
        data.iter().map(|x| x * i16::MAX as f64).collect()
    };

    // plot waveform, if necessary, handle errors
    if let Some(path) = plot_waveform {
        plot_numbers(&path, data)?
    }

    // write .wav file, handle errors
    for sample in samples.into_iter() {
        writer.write_sample(sample as i16)?
    }

    Ok(())
}

/// reads the sample data into a Vector, the channel data is interleaved.
pub fn read_from_wav(file_path: &str) -> Vec<f64> {
    let mut reader = hound::WavReader::open(file_path).unwrap();
    reader.samples::<i32>().map(| s | s.unwrap() as f64).collect()
}

/// returns a WavSpec struct with fields: channels, sample_rate, bits_per_sample, sample_format
pub fn get_wav_specs(file_path: &str) -> WavSpec {
    hound::WavReader::open(file_path).unwrap().spec()
}