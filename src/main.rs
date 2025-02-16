use std::path::Path;
use image_dft::*;
use std::process;
use clap::Parser;

fn main() {
    // parse command line arguments
    let cli = Cli::parse();

    // make sure an argument for filepath was supplied
    let file_path: &str = &cli.input.clone().unwrap_or_else(|| {
        println!("No input file!");
        process::exit(1);
    });
    let file_ext = Path::new(file_path).extension().and_then(|s| s.to_str()).unwrap();

    // file extension determines, what function is called
    if ["jpg", "png"].contains(&file_ext) {
        image_to_waves(file_path, &cli)
    } else if ["wav"].contains(&file_ext) {
        // transform_sound_to_image(file_path, table_len)
        sound_to_img_sequence(file_path, &cli)
    } else {
        Err(format!("Sorry, this file format is not supported at the moment: {}", file_ext).into())
    }
        .unwrap_or_else(|err| {
            println!("{}", err);
            process::exit(1);
        });
}