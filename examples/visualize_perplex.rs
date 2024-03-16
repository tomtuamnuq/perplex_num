use perplex_num::Perplex;
use plotters::prelude::*;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let root = BitMapBackend::new("./examples/perplex_numbers.jpg", (640, 480)).into_drawing_area();
    root.fill(&WHITE)?;
    let (x_min, x_max, y_min, y_max) = (-1, 4, -1, 2);

    let mut chart = ChartBuilder::on(&root)
        .caption(
            "Perplex Number Visualization",
            ("sans-serif", 50).into_font(),
        )
        .margin(5)
        .x_label_area_size(30)
        .y_label_area_size(30)
        .build_cartesian_2d(x_min as f64..x_max as f64, y_min as f64..y_max as f64)?;

    chart.configure_mesh().draw()?;

    let z = Perplex { t: 1.0, x: 0.5 };
    let font = ("sans-serif", 15).into_font();
    let label_offset = (10, 0);
    let numbers = [("z", z), ("conj", z.conj()), ("h", Perplex::h())];
    let functions = [
        ("inv", z.try_inverse().expect("z is invertible")),
        ("sinus", z.sin()),
        ("sinh", z.sinh()),
        ("exp", z.exp()),
        ("ln", z.ln().expect("The natural logarithm of z exists!")),
    ];
    // Draw a diagonal line where x = y
    chart
        .draw_series(LineSeries::new(
            (y_min..=y_max).map(|x| (x as f64, x as f64)),
            &BLACK.mix(0.75),
        ))?
        .label("Light-like numbers t=x")
        .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], BLACK));
    chart.draw_series(LineSeries::new(
        (y_min..y_max).map(|x| (x as f64, -x as f64)),
        &BLACK.mix(0.75),
    ))?;
    // Draw the hyperbola in the right section
    let d = z.squared_distance();
    let mut hyperbola = Vec::new();
    (0..=1000)
        .map(|i| d + i as f64 * (x_max as f64 - d) / 1000.0)
        .for_each(|x| {
            let y = (x * x - d).sqrt();
            if y < y_max as f64 {
                hyperbola.push((x, y));
            }
            let y = -y;
            if y > y_min as f64 {
                hyperbola.push((x, y));
            }
        });
    chart
        .draw_series(
            hyperbola
                .into_iter()
                .map(|coord| Circle::new(coord, 1, RED.mix(0.5))),
        )?
        .label(format!("Hyperbola defined by t²-x²={:.2}", d))
        .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], RED));
    // Plot the original Perplex numbers
    chart.draw_series(numbers.iter().map(|(label, z)| {
        let coord = (z.t, z.x);
        EmptyElement::at(coord)
            + Circle::new((0, 0), 4, RED.filled())
            + Text::new(
                format!("{}({:.2}, {:.2})", label, coord.0, coord.1),
                label_offset,
                font.clone(),
            )
    }))?;
    // Plot the results of the functions
    chart.draw_series(functions.into_iter().map(|(label, z)| {
        let coord = (z.t, z.x);
        EmptyElement::at(coord)
            + Circle::new((0, 0), 5, BLUE.filled())
            + Text::new(
                format!("{}({:.2}, {:.2})", label, coord.0, coord.1),
                label_offset,
                font.clone(),
            )
    }))?;
    chart.configure_series_labels().draw()?;
    Ok(())
}
