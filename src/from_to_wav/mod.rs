use hound::WavSpec;
use crate::utilities::rescale;
use crate::plot::plot_numbers;

// TODO normalize should use the max of absolute numbers
/// Write numbers from a vector into a .wav file. Numbers should be between -1.0 and 1.0,
/// else, when normalize is true, they are rescaled to fit that range.
pub fn write_to_wav (name: &str, data: &[f64], normalize: bool, plot_waveform: bool) -> Result<(), Box<dyn std::error::Error>>  {

    let spec = WavSpec {
        channels: 1,
        sample_rate: 48000,
        bits_per_sample: 16,
        sample_format: hound::SampleFormat::Int,
    };

    let mut writer = hound::WavWriter::create(name, spec).unwrap();

    // normalize if necessary
    let samples: Vec<f64> =  if normalize {
        let max = data.iter().fold(0.0, |a: f64, &b| a.max(b.abs()));
        let min = data.iter().fold(f64::INFINITY, |a: f64, &b| a.min(b.abs()));

        data.iter().map(|x| rescale(*x, min, max, -(i16::MAX as f64), i16::MAX as f64)).collect()
    } else {
        data.iter().map(|x| x * i16::MAX as f64).collect()
    };

    // plot waveform, if necessary, handle errors
    if plot_waveform {
        plot_numbers("waveform.png", data)?
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