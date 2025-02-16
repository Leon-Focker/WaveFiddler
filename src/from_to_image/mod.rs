use image::{ImageBuffer, RgbImage};
use std::error::Error;

pub fn write_to_image(name: &str, pixel_data: &[u8], width: usize, height: usize) -> Result<(), Box<dyn Error>> {
    // Create an RgbImage from the pixel data
    let image: RgbImage = ImageBuffer::from_vec(width as u32, height as u32, pixel_data.to_vec())
        .ok_or_else(|| "Buffer size does not match image dimensions".to_string())?; // Proper error conversion

    // Save the image to a file
    image.save(name).map_err(Into::into) // Converts ImageError into Box<dyn Error>
}
