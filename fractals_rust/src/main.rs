extern crate image;
extern crate num_complex;
extern crate num_integer;
extern crate rgsl;
extern crate rayon;

use rayon::prelude::*;

struct Image{
    imgx: u32,
    imgy: u32,
    //imgbuff: image::ImageBuffer,
}

// On veut pouvoir definir des nouvelles fractales facilement
struct Fractale <T>
    where T: Fn(num_complex::Complex<f64>, num_complex::Complex<f64>
                ) -> num_complex::Complex<f64> + Sync
{
    rangex: (f64, f64),
    rangey: (f64, f64),
    equation: T,
    max_iter: u32,
}

impl <T> Fractale <T>
    where T: Fn(num_complex::Complex<f64>, num_complex::Complex<f64>
                ) -> num_complex::Complex<f64> + Sync
{
    fn recenter (&mut self, new_ctr :(f64, f64)) {
        let centrex = (self.rangex.0 + self.rangex.1) / 2.0;
        let centrey = (self.rangey.0 + self.rangey.1) / 2.0;
        self.rangex.0 += (new_ctr.0 - centrex).abs();
        self.rangex.1 += (new_ctr.0 - centrex).abs();
        self.rangey.0 += (new_ctr.1 - centrey).abs();
        self.rangey.1 += (new_ctr.1 - centrey).abs();
    }

    fn zoom (&mut self, zoom: f64) {
        let spanx = (self.rangex.1 - self.rangex.0).abs();
        let spany = (self.rangey.1 - self.rangey.0).abs();
        
        let centrex = (self.rangex.0 + self.rangex.1) / 2.0;
        let centrey = (self.rangey.0 + self.rangey.1) / 2.0;
        self.rangex.0 = centrex - (spanx/2.0) / zoom;
        self.rangex.1 = centrex + (spanx/2.0) / zoom;
        self.rangey.0 = centrey - (spany/2.0) / zoom;
        self.rangey.1 = centrey + (spany/2.0) / zoom;
    }
}

fn make_image<T>(frac: &Fractale<T>, min_size: u32) -> Image
    where T: Fn(num_complex::Complex<f64>, num_complex::Complex<f64>
                ) -> num_complex::Complex<f64> + Sync
    {
        // trouver max_size en fonction de rangex et rangey
        let spanx = (frac.rangex.1 - frac.rangex.0).abs();
        let spany = (frac.rangey.1 - frac.rangey.0).abs();
        let span_max = spanx.max(spany);
        let span_min = spanx.min(spany);
        let max_size = ((min_size as f64) * span_max / span_min) as u32;

        // imgx = max si spanx > spany
        Image {
            imgx: if spanx > spany {max_size} else {min_size},
            imgy: if spanx > spany {min_size} else {max_size},
        }
    }

fn true_hist <T> (img: Image, frac: Fractale<T>) 
    // Eviter de copier / coller T tout le temps!
    where T: Fn(num_complex::Complex<f64>, num_complex::Complex<f64>
                ) -> num_complex::Complex<f64> + Sync
    {
    // Create a new ImgBuf with width: imgx and height: imgy
    let mut imgbuf: image::RgbImage = image::ImageBuffer::new(img.imgx, img.imgy);
    // scaling factors
    let scalex = (frac.rangex.1 - frac.rangex.0) / img.imgx as f64;
    let scaley = (frac.rangey.1 - frac.rangey.0) / img.imgy as f64;
    // iterate over pixels
    let pixels: Vec<f64> = (0..img.imgx * img.imgy)
        .into_par_iter()
        .map( |j| {
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
            let mut mu:f64 = -1.0;
            let mut i = 0;
            while i < frac.max_iter-1 {
                z = (frac.equation)(z,c);
                i += 1;
                if z.norm() > 400.0 {
                    mu = (i + 1) as f64 - z.norm().log(10.0).log(2.0);
                    break;
                }
            }
            mu
        }).collect();
    let mut histogram = pixels.clone();
    // compute histogram
    histogram.sort_by(|a, b| a.partial_cmp(b).unwrap());
    // find index of interior pixels
    let mut interior:u32 = 0;
    for i in 0..img.imgx * img.imgy {
        if histogram[i as usize] != -1.0 {
            interior = i;
            break;
        }
    }
    let exterior:f64 = 1.0 / ((img.imgx * img.imgy - interior) as f64);
    // lookup histogram
    
    // mettre un cache ?
    // on veut un tableau qui retourne le rang dans l'histo de chaque mu 
    // du tableau pixels
    // attention ! histogram n'est pas vraiment un histogramme
    // parcourir histogram et à chaque changement de valeur 
    // enregister un dic {val; idx}. ensuite lookup par val
    //let ix = histogram.iter().position(|&r| r == mu as f64).unwrap() as u32;
    let grey_vec: Vec<f64> = (0..img.imgx * img.imgy)
        .into_par_iter()
        .map(|i| {
            let mu = pixels[i as usize] as f64;
            let ix = rgsl::interpolation::bsearch(
                &histogram[..], mu, 0, (img.imgx * img.imgy) as usize) as u32;
            let mut grey = 0.0;
            if interior <= ix {
                //grey =((ix - interior) as f64 * exterior).powf(3.0);
                let x = ((ix - interior) as f64 * exterior);
                grey = 1.0 - (1.0 - x).powf(1.0/1.5);
            }
            grey
        }).collect();

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
    
    let mandel = Fractale {
        rangex: (-2.5, 1.0),
        rangey: (-1.3, 1.3),
        equation: |z,c| z * z + c,
        max_iter: 5000,
    };
    let tricorn = Fractale {
        rangex: (-2.0, 2.0),
        rangey: (-2.0, 2.0),
        equation: |z,c| z.conj() * z.conj() + c,
        max_iter: 30,
    };

    let mut frac = mandel;

    // zoom et recadrage
    let new_ctr = (-0.787, 0.25);
    let zoom = 100.0;
    frac.recenter(new_ctr);
    frac.zoom(zoom);
    // generer la taille d'image
    let min_size = 1000;
    let img = make_image(&frac, min_size);
    
    // compute !
    true_hist(img, frac);
}













/*
fn palette (x: f32, r:f32, g:f32, b:f32) -> [u8; 3] {
    let c = 1.0 - x;
    [(c.powf(r)*255.0) as u8,  (c.powf(g)*255.0) as u8, (c.powf(b)*255.0) as u8]
}

fn classic <T> (img: Image, frac: Fractale<T>)
    where T: Fn(num_complex::Complex<f32>, num_complex::Complex<f32>
                ) -> num_complex::Complex<f32>
{
    // Create a new ImgBuf with width: imgx and height: imgy
    let mut imgbuf: image::RgbImage = image::ImageBuffer::new(img.imgx, img.imgy);

    /* 1st pass
     * Calculate iteration counts for each pixel
     */
    let scalex = (frac.rangex.1 - frac.rangex.0) / img.imgx as f32;
    let scaley = (frac.rangey.1 - frac.rangey.0) / img.imgy as f32;

    for x in 0..img.imgx {
        for y in 0.. img.imgy {
            let cx = x as f32 * scalex + frac.rangex.0;
            let cy = y as f32 * scaley + frac.rangey.0;
            let c = num_complex::Complex::new(cx, cy);
            let mut z = num_complex::Complex::new(0.0, 0.0);

            let mut i = 0;
            while i < frac.max_iter && z.norm() <= 2.0 {
                z = (frac.equation)(z,c);
                i += 1;
            }
            // NAN quand z.norm() < 1
            let up = (z.norm().log(10.0).log(10.0))
                /((frac.max_iter as f32).log(10.0));
            let s = i as f32 - up / (2.0f32.log(10.0));
            let pixel = imgbuf.get_pixel_mut(x,y);
            let image::Rgb(data) = *pixel;
            *pixel = image::Rgb(palette(s / frac.max_iter as f32, 2.0, 1.5, 1.0));
        }
    }
    imgbuf.save("fractal.png").unwrap();
}

fn histogram <T> (img: Image, frac: Fractale<T>) 
    // Eviter de copier / coller T tout le temps!
    where T: Fn(num_complex::Complex<f32>, num_complex::Complex<f32>
                ) -> num_complex::Complex<f32>
    {
    // Create a new ImgBuf with width: imgx and height: imgy
    let mut imgbuf: image::RgbImage = image::ImageBuffer::new(img.imgx, img.imgy);

    /* 1st pass
     * Calculate iteration counts for each pixel
     */
    let scalex = (frac.rangex.1 - frac.rangex.0) / img.imgx as f32;
    let scaley = (frac.rangey.1 - frac.rangey.0) / img.imgy as f32;
    let mut itcount = vec![vec![0u32; img.imgy as usize]; img.imgx as usize];

    for x in 0..img.imgx {
        for y in 0.. img.imgy {
            let cx = x as f32 * scalex + frac.rangex.0;
            let cy = y as f32 * scaley + frac.rangey.0;
            let c = num_complex::Complex::new(cx, cy);
            let mut z = num_complex::Complex::new(0.0, 0.0);

            let mut i = 0;
            while i < frac.max_iter-1 && z.norm() <= 2.0 {
                z = (frac.equation)(z,c);
                i += 1;
            }
            itcount[x as usize][y as usize] = i;
        }
    }
    /* 2nd pass
     * Create an array histo of size n = max iteration count
     * Iterate over itcount and build the histogram
     */
    let mut histo = vec![0u32; frac.max_iter as usize];

    for x in 0..img.imgx {
        for y in 0.. img.imgy {
            let i = itcount[x as usize][y as usize];
            histo[i as usize] += 1;
        }
    }
    /* 3rd pass
     * Count the number of iterations */
    let mut total: u32 = 0;
    for i in 0..frac.max_iter {
        total += histo[i as usize];
    }
    /* 4th pass
     * All values in itcounts are indexed each itcount i is normalized
     */
    let mut hue = vec![vec![0f32; img.imgy as usize]; img.imgx as usize];
    for x in 0..img.imgx {
        for y in 0..img.imgy {
            let iteration = itcount[x as usize][y as usize];
            for i in 0..iteration {
                hue[x as usize][y as usize] += histo[i as usize] as f32 / total as f32;
            }
        }
    }
    for (x, y, pixel) in imgbuf.enumerate_pixels_mut() {
        /* Print image */
        let image::Rgb(data) = *pixel;
        *pixel = image::Rgb(palette(hue[x as usize][y as usize], 1.0, 1.0, 1.0));
    }
    // Save the image as “fractal.png”, the format is deduced from the path
    imgbuf.save("fractal.png").unwrap();
    println!("{:?}", histo);
}
*/
