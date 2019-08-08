# fractals_rust
A simple parallel fractal render in Rust, writen as a practice exercise.

For now this supports :
- Zooming and reframing
- Custom resolution
- The Mandelbrot and Tricorn fractals

It uses the histogram equalization algorithm.
The algorithm is roughly the following :

```
1. Make two vectors. One for the saving the iteration counts, and one to compute a "histogram" (actually a vector of sorted pixels).
2. Generate the image using the smooth coloring algorithm[1].
3. Copy the generated image into the histogram vector.
4. Sort the histogram
5. Lookup in the histogram the escape iteration count corresponding to the inside of the fractal, ie the divergence point
6. Color the pixels according to their distance from the divergence point
```

Tip : when doing big zooms, increase the iteration count to have sharper images.

To run : `cargo run --release`

Improvement idea : instead of a vector for the histogram, use a Hashmap storing (iteration count, number of pixels) (ie an actual histogram). We may have to do some bining and thus reduce the color resolution though. However that would put the lookup time from O(log(n)) (binary search) to O(1) while reducing memory usage (depending on the size of the bins)

### References
[1] https://en.wikipedia.org/wiki/Mandelbrot_set#Continuous_(smooth)_coloring
