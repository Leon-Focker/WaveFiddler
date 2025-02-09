
pub fn rescale<T: Copy + std::ops::Add<Output = T> + std::ops::Sub<Output = T> + std::ops::Div<Output = T>+ std::ops::Mul<Output = T>>
(value: T, old_min: T, old_max: T, new_min: T, new_max: T) -> T {
    (((value - old_min) / (old_max - old_min)) * (new_max - new_min)) + new_min
}

// copied from fft_2d because it is not a public function:
pub fn transpose<T: Copy + Default>(width: usize, height: usize, matrix: &[T]) -> Vec<T> {
    let mut ind = 0;
    let mut ind_tr;
    let mut transposed = vec![T::default(); matrix.len()];
    for row in 0..height {
        ind_tr = row;
        for _ in 0..width {
            transposed[ind_tr] = matrix[ind];
            ind += 1;
            ind_tr += height;
        }
    }
    transposed
}

fn partial_image_32<T: Copy>(image_data: &Vec<T>, width: usize, height: usize, (start_x, start_y): (usize, usize)) -> Result<Vec<T>, ()> {
    if start_x + 32 >= width || start_y + 32 >= height || width * height > image_data.len() {
       return Err(())
    }

    let mut result = Vec::with_capacity(32 * 32);

    // go through 32 rows
    for y in start_y..(start_y + 32) {
        let x = start_x + (y * width as usize);
        for idx in x..(x + 32) {
            result.push(image_data[idx])
        }
    };

    Ok(result)
}

pub fn split_image<T: Copy>(image_data: &Vec<T>, width: usize, height: usize) -> Vec<Vec<T>> {
    let mut result: Vec<Vec<T>> = Vec::new();
    let per_row = width / 32;
    let per_column = height / 32;

    for x in 0..per_row {
        for y in 0..per_column {
            result.push(partial_image_32(image_data, width, height, (x * 32, y * 32)).unwrap())
        }
    }

    result
}