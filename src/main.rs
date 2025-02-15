use std::path::Path;
use image_dft::*;
use clap::Parser;

/// Command Line Arguments
#[derive(Parser)]
#[command(version, about, long_about = None)]
struct Cli {
    /// Path to the input file
    input: Option<String>,

    /// Window Size the fft, defaults to 512
    #[arg(short, long, default_value_t = 512)]
    fft_size: u32,
}

fn main() -> Result<(), Box<dyn std::error::Error>>{
    // parse command line arguments
    let cli = Cli::parse();

    if cli.input.is_none() {
        panic!("no input file!")
    }

    let file_path: &str = &cli.input.unwrap();
    let table_len: usize = cli.fft_size.try_into().unwrap();
    let file_ext = Path::new(file_path).extension().and_then(|s| s.to_str()).unwrap();

    // file extension determines, what function is called
    if ["jpg", "png"].contains(&file_ext) {
        transform_image_to_sound(file_path, table_len)
    } else if ["wav"].contains(&file_ext) {
        // transform_sound_to_image(file_path, table_len)
        sound_to_img_sequence(file_path)
    } else {
        println!("sorry, this file format is not supported at the moment: {}", file_ext);
        Ok(())
    }
}