extern crate image;
extern crate num_complex;
extern crate num_integer;
extern crate rayon;

use rayon::prelude::*;

/// Represents the output image resolution
struct Image {
    imgx: u32,
    imgy: u32,
}

/// We define a Fractal struct for defining new
/// fractal equations easily.
/// - rangex : the range of the fractal in the real axis
///- rangey : the range of the fractal on the complex axis
/// - equation : a closure defining the fractal equation
/// - max_iter : the maximum amount of iterations allowed per pixel
struct Fractal<T>
where
    T: Fn(num_complex::Complex<f64>, num_complex::Complex<f64>) -> num_complex::Complex<f64> + Sync,
{
    rangex: (f64, f64),
    rangey: (f64, f64),
    equation: T,
    max_iter: u32,
}

/// Helper functions for recentering and zooming on the fractal.
///  Works by adding offsets or scaling to the ranges
impl<T> Fractal<T>
where
    T: Fn(num_complex::Complex<f64>, num_complex::Complex<f64>) -> num_complex::Complex<f64> + Sync,
{
    fn recenter(&mut self, new_ctr: (f64, f64)) {
        let centrex = (self.rangex.0 + self.rangex.1) / 2.0;
        let centrey = (self.rangey.0 + self.rangey.1) / 2.0;
        self.rangex.0 += (new_ctr.0 - centrex).abs();
        self.rangex.1 += (new_ctr.0 - centrex).abs();
        self.rangey.0 += (new_ctr.1 - centrey).abs();
        self.rangey.1 += (new_ctr.1 - centrey).abs();
    }

    fn zoom(&mut self, zoom: f64) {
        let spanx = (self.rangex.1 - self.rangex.0).abs();
        let spany = (self.rangey.1 - self.rangey.0).abs();

        let centrex = (self.rangex.0 + self.rangex.1) / 2.0;
        let centrey = (self.rangey.0 + self.rangey.1) / 2.0;
        self.rangex.0 = centrex - (spanx / 2.0) / zoom;
        self.rangex.1 = centrex + (spanx / 2.0) / zoom;
        self.rangey.0 = centrey - (spany / 2.0) / zoom;
        self.rangey.1 = centrey + (spany / 2.0) / zoom;
    }
}

/// For a given fractal and minimum size in pixels, returns an image
/// struct of the correct aspect ratio
fn make_image<T>(frac: &Fractal<T>, min_size: u32) -> Image
where
    T: Fn(num_complex::Complex<f64>, num_complex::Complex<f64>) -> num_complex::Complex<f64> + Sync,
{
    // find max_size according to rangex and rangey
    let spanx = (frac.rangex.1 - frac.rangex.0).abs();
    let spany = (frac.rangey.1 - frac.rangey.0).abs();
    let span_max = spanx.max(spany);
    let span_min = spanx.min(spany);
    let max_size = ((min_size as f64) * span_max / span_min) as u32;

    // imgx = max if spanx > spany
    Image {
        imgx: if spanx > spany { max_size } else { min_size },
        imgy: if spanx > spany { min_size } else { max_size },
    }
}

/// Generates the fractal
/// Uses an histogram equalization technique involving
/// sorting the image pixels (with quicksort) and then looking
/// through it (with binary search).
/// Complexity : (n^2 log(n)), n = side of the image */
fn true_hist<T>(img: Image, frac: Fractal<T>)
// Is there any way to avoid copy / pasting this all the time ?
where
    T: Fn(num_complex::Complex<f64>, num_complex::Complex<f64>) -> num_complex::Complex<f64> + Sync,
{
    // Create a new ImgBuf with width: imgx and height: imgy
    let mut imgbuf: image::RgbImage = image::ImageBuffer::new(img.imgx, img.imgy);
    // scaling factors
    let scalex = (frac.rangex.1 - frac.rangex.0) / img.imgx as f64;
    let scaley = (frac.rangey.1 - frac.rangey.0) / img.imgy as f64;
    // iterate over pixels
    let pixels: Vec<f64> = (0..img.imgx * img.imgy)
        .into_par_iter()
        .map(|j| {
            // scale the pixels
            let xy = num_integer::div_rem(j, img.imgy);
            let x = xy.0;
            let y = xy.1;
            let cx = x as f64 * scalex + frac.rangex.0;
            let cy = y as f64 * scaley + frac.rangey.0;
            // define the fractal
            let c = num_complex::Complex::new(cx, cy);
            let mut z = num_complex::Complex::new(0.0, 0.0);
            // smooth counter
            let mut mu: f64 = -1.0;
            let mut i = 0;
            while i < frac.max_iter - 1 {
                z = (frac.equation)(z, c);
                i += 1;
                if z.norm() > 400.0 {
                    mu = (i + 1) as f64 - z.norm().log(10.0).log(2.0);
                    break;
                }
            }
            mu
        })
        .collect();
    let mut histogram = pixels.clone();
    // compute histogram
    histogram.sort_by(|a, b| a.partial_cmp(b).unwrap());
    // find index of interior pixels
    let mut interior: u32 = 0;
    for i in 0..img.imgx * img.imgy {
        if histogram[i as usize] != -1.0 {
            interior = i;
            break;
        }
    }
    let exterior: f64 = 1.0 / ((img.imgx * img.imgy - interior) as f64);
    // lookup histogram
    // we want to look for each pixel how far away from the interior it is
    // from an iteration standpoint (not spatial)
    let grey_vec: Vec<f64> = (0..img.imgx * img.imgy)
        .into_par_iter()
        .map(|i| {
            let mu = pixels[i as usize] as f64;
            let ix = match histogram.binary_search_by(|a| a.partial_cmp(&mu).unwrap()) {
                Ok(v) => v as u32,
                Err(e) => std::cmp::min(e, (img.imgx * img.imgy) as usize - 1) as u32,
            };
            let mut grey = 0.0;
            if interior <= ix {
                let x = ((ix - interior) as f64 * exterior);
                grey = 1.0 - (1.0 - x).powf(1.0 / 1.5);
            }
            grey
        })
        .collect();

    for i in 0..img.imgx * img.imgy {
        let grey = grey_vec[i as usize];
        // write image
        let tmp = num_integer::div_rem(i, img.imgy);
        let pixel = imgbuf.get_pixel_mut(tmp.0, tmp.1);
        let image::Rgb(data) = *pixel;
        let g = (255.0 * grey) as u8;
        *pixel = image::Rgb([g, g, g]);
    }
    // Save the image as “fractal.png”, the format is deduced from the path
    imgbuf.save("fractal.png").unwrap();
}

fn main() {
    let mandel = Fractal {
        rangex: (-2.5, 1.0),
        rangey: (-1.3, 1.3),
        equation: |z, c| z * z + c,
        max_iter: 5000,
    };
    let tricorn = Fractal {
        rangex: (-2.0, 2.0),
        rangey: (-2.0, 2.0),
        equation: |z, c| z.conj() * z.conj() + c,
        max_iter: 30,
    };

    let mut frac = mandel;

    // recenter and zoom
    // the "normal" plot for Mandelbrot is new_ctr = (-0.75, 0) and zoom = 1.0
    let new_ctr = (-0.787, 0.25);
    let zoom = 100.0;
    frac.recenter(new_ctr);
    frac.zoom(zoom);
    // make image canvas
    let min_size = 1000;
    let img = make_image(&frac, min_size);

    // compute !
    true_hist(img, frac);
}
