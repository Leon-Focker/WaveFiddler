use rustfft::num_complex::Complex;
use plotters::prelude::*;

pub fn _plot_complex_numbers(input: &[Complex<f64>]) -> Result<(), Box<dyn std::error::Error>> {
    let x_min = input.iter().fold(f64::INFINITY, |a, &b| a.min(b.re));
    let x_max = input.iter().fold(f64::NEG_INFINITY, |a, &b| a.max(b.re));
    let y_min = input.iter().fold(f64::INFINITY, |a, &b| a.min(b.im));
    let y_max = input.iter().fold(f64::NEG_INFINITY, |a, &b| a.max(b.im));

    // Create a new chart
    let root = BitMapBackend::new("complex_plot.png", (640, 480)).into_drawing_area();
    root.fill(&WHITE)?;

    // Set up the chart
    let mut chart = ChartBuilder::on(&root)
        .caption("Simple Vector Plot", ("Arial", 40))
        .build_cartesian_2d(x_min..x_max, y_min.round()..y_max.round())?;

    // start with origin
    chart.draw_series(LineSeries::new(vec![(0.0, 0.0)].into_iter(), &BLUE).point_size(5))?;

    // Draw the x-y scatter plot
    chart.draw_series(LineSeries::new(
        input.iter().map(|complex| (complex.re, complex.im)),
        &RED,
    ))?;

    Ok(())
}

pub fn plot_numbers(name: &str, input: &[f64]) -> Result<(), Box<dyn std::error::Error>> {
    let y_min = input.iter().fold(f64::INFINITY, |a, &b| a.min(b));
    let y_max = input.iter().fold(f64::NEG_INFINITY, |a, &b| a.max(b));

    // Create a new chart
    let root = BitMapBackend::new(name, (640, 480)).into_drawing_area();
    root.fill(&WHITE)?;

    // Set up the chart
    let mut chart = ChartBuilder::on(&root)
        .caption(name, ("Arial", 40))
        .build_cartesian_2d(0..(input.len() as i32), y_min.round()..y_max.round())?;

    // Draw the x-y scatter plot
    chart.draw_series(LineSeries::new(
        input.iter().enumerate().map(|(x, y)| (x as i32, *y)),
        &RED,
    ))?;

    Ok(())
}