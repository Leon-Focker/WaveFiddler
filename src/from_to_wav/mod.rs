use crate::utilities::rescale;
use crate::plot::plot_numbers;

pub fn write_to_wav (name: &str, data: &[f64], normalize: bool) -> Result<(), Box<dyn std::error::Error>>  {

    let spec = hound::WavSpec {
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

    let _ = plot_numbers("wavetable.png", &samples);

    for sample in samples.into_iter() {
        writer.write_sample(sample as i16).unwrap();
    }

    Ok(())
}

/// reads the sample data into a Vector
pub fn read_from_wav(file_path: &str) -> Vec<f64> {
    let mut reader = hound::WavReader::open(file_path).unwrap();
    reader.samples::<i32>().map(| s | s.unwrap() as f64).collect()
}