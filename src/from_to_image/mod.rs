use image::{ImageBuffer, RgbImage};

pub fn write_to_image(name: &str, pixel_data: &[u8], width: usize, height: usize) -> Result<(), Box<dyn std::error::Error>> {

    // Create an RgbImage from the pixel data
    let image: RgbImage = ImageBuffer::from_vec(width as u32, height as u32, pixel_data.to_vec())
        .expect("Buffer size does not match image dimensions");

    // Save the image to a file
    image.save(name).expect("Failed to save image");

    Ok(())
}