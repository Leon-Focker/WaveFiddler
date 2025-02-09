use std::i16;
use hound;
use crate::utilities::rescale;
use crate::plot::plot_numbers;

pub fn write_to_wav (name: &str, data: &Vec<f64>) -> Result<(), Box<dyn std::error::Error>>  {

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