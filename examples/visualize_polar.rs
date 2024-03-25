use perplex_num::{HyperbolicPolar, Perplex};
use plotters::{
    prelude::*,
    style::full_palette::{LIGHTBLUE, LIGHTGREEN},
};
use std::iter;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // Create a root drawing area and split it into two
    let root = BitMapBackend::new("./examples/polar.jpg", (1200, 600)).into_drawing_area();
    root.fill(&WHITE)?;
    let (left, right) = root.split_horizontally(600);

    // Define the range for t and x
    let (t_min, t_max, x_min, x_max) = (0, 3, -2, 2);
    // Define the range for rho theta
    let (rho_min, rho_max, theta_min, theta_max) = (0, 2, -2, 2);
    // Create the left chart builder for the standard form
    let mut left_chart = ChartBuilder::on(&left)
        .caption(
            "Perplex Numbers in Standard Form",
            ("sans-serif", 30).into_font(),
        )
        .margin(30)
        .x_label_area_size(100)
        .y_label_area_size(100)
        .build_cartesian_2d(t_min as f64..t_max as f64, x_min as f64..x_max as f64)?;
    // Create the right chart builder for the polar form
    let mut right_chart = ChartBuilder::on(&right)
        .caption(
            "Perplex Numbers in Hyperbolic Polar Form",
            ("sans-serif", 30).into_font(),
        )
        .margin(30)
        .x_label_area_size(100)
        .y_label_area_size(100)
        .build_cartesian_2d(
            rho_min as f64..rho_max as f64,
            theta_min as f64..theta_max as f64,
        )?;

    let font = ("sans-serif", 20).into_font();
    left_chart
        .configure_mesh()
        .max_light_lines(2)
        .label_style(font.clone())
        .x_labels(3)
        .y_labels(3)
        .x_label_formatter(&|t| format!("t={:.0}", t))
        .y_label_formatter(&|x| format!("x={:.0}", x))
        .draw()?;

    right_chart
        .configure_mesh()
        .max_light_lines(2)
        .label_style(font.clone())
        .x_labels(5)
        .y_labels(3)
        .x_label_formatter(&|rho| format!("ρ={:.1}", rho))
        .y_label_formatter(&|theta| format!("θ={:.0}", theta))
        .draw()?;
    let down_alpha = 0.1;
    let z = Perplex { t: 1.0, x: 0.5 };
    let d = z.squared_distance();
    let mut hyperbola_ru = Vec::new();
    let mut hyperbola_rd = Vec::new();
    (0..=100_000)
        .map(|i| d + i as f64 * (t_max as f64 - d) / 100_000.0)
        .for_each(|t| {
            let x = (t * t - d).sqrt();
            hyperbola_ru.push(Perplex::new(t, x));
            hyperbola_rd.push(Perplex::new(t, -x));
        });

    // Highlight the Right sector
    left_chart.draw_series(LineSeries::new(
        (t_min..=t_max).map(|x| (x as f64, x as f64)),
        &BLACK,
    ))?;
    left_chart.draw_series(LineSeries::new(
        (t_min..t_max).map(|x| (x as f64, -x as f64)),
        &BLACK,
    ))?;
    let bound_filter_perplex = |z: &Perplex<f64>| {
        t_min as f64 <= z.t && z.t <= t_max as f64 && x_min as f64 <= z.x && z.x <= x_max as f64
    };
    left_chart.draw_series(iter::once(Text::new(
        format!("t²-x²={:.2}", d),
        (3.0 * d, 0.1),
        font.clone().color(&BLUE),
    )))?;
    let _bound_filter_polar = |p: &HyperbolicPolar<f64>| {
        rho_min as f64 <= p.rho
            && p.rho <= rho_max as f64
            && theta_min as f64 <= p.theta
            && p.theta <= theta_max as f64
    };
    let perplex_coords = |z: &Perplex<f64>| (z.t, z.x);
    let polar_coords = |p: &HyperbolicPolar<f64>| (p.rho, p.theta);
    let (mut result_perplex, mut result_polar) = (Vec::new(), Vec::new());
    for (i, color) in [
        (-3, &LIGHTGREEN),
        (-2, &MAGENTA),
        (-1, &LIGHTBLUE),
        (1, &BLUE),
        (2, &RED),
        (3, &GREEN),
    ] {
        for (hyperbola, alpha) in [(&hyperbola_rd, down_alpha), (&hyperbola_ru, 1.0)] {
            result_perplex.clear();
            result_polar.clear();
            for z in hyperbola.iter() {
                if let Some(z_pow) = z.powi(i) {
                    if bound_filter_perplex(&z_pow) {
                        let polar = z_pow.polar();
                        result_perplex.push(z_pow);
                        result_polar.push(polar);
                    }
                }
            }
            left_chart.draw_series(LineSeries::new(
                result_perplex.iter().map(perplex_coords),
                color.mix(alpha),
            ))?;
            right_chart.draw_series(LineSeries::new(
                result_polar.iter().map(polar_coords),
                color.mix(alpha),
            ))?;
        }

        let z_middle = result_perplex[result_perplex.len() / 2];
        let label_pos_perplex = perplex_coords(&z_middle);
        let polar_middle = result_polar[result_polar.len() / 2];
        let label_pos_polar = polar_coords(&polar_middle);
        let _left_draw_res = left_chart.draw_series(iter::once(Text::new(
            format!("{i}"),
            (label_pos_perplex.0 + 0.1, 1.0),
            font.clone().color(color),
        )))?;
        let right_draw_res = right_chart.draw_series(iter::once(Text::new(
            format!("{i}"),
            (label_pos_polar.0 + 0.01, 1.0),
            font.clone().color(color),
        )))?;
        let legend = match i {
            1 => format!("Base Hyperbola H at t²-x²={:.2}", d),
            2 => String::from("H to the power of 2 etc."),
            _ => String::from(""),
        };
        if !legend.is_empty() {
            right_draw_res
                .label(legend)
                .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], color.to_owned()));
        }
    }
    // Configure and draw the series labels
    for chart in [&mut left_chart, &mut right_chart].iter_mut() {
        chart
            .configure_series_labels()
            .position(SeriesLabelPosition::UpperRight)
            .label_font(font.clone())
            .draw()?;
    }
    Ok(())
}
