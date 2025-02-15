use std::path::Path;
use std::env;
use image_dft::*;

fn main() -> Result<(), Box<dyn std::error::Error>>{

    let args: Vec<String> = env::args().collect();

    // important parameters:
    // TODO error when not supplied
    let file_path: &str = &args[1];
    let table_len: usize = args[2].parse::<usize>().unwrap_or(512);
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