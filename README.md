# PALETTE
**PALETTE**: *Perceived Average of Light in Each Temporal Tile Extracted*

PALETTE extracts representative color palettes from visual data using k-means clustering in perceptually uniform color spaces.

## Current Status ##
- **SoA Migration Complete**: 'struct Color' now stores contiguous channel arrays ('alpha', 'beta', 'gamma') with 32-byte alignment for direct vector SIMD loads.

- **Color Space Conversions**: Bidirectional conversions between linear sRGB, Oklab, and CIELAB with SSE and AVX2 routines where available.

- **Distance Metric**: Euclidean distance, DIE94, and CIEDE2000 implementations.

- **KMeans Clustering**: Full k-means++ centroid initialization and clustering pipeline on SoA color batches.

- **Frame & Image Processing**:
  - Full FFmpeg video frame decoding and YUV $/to$ RGB scaling via `libswscale` into linear sRGB SoA.
  - Standalone PPM (P3/P6) and BMP image loader for dependency-free palette extraction from images.
  - Palette swatch image exporter (`.ppm`).

- **Cross-Platform**: Build natively on macOS (Apple Silicon / Intel), Linux, and Windows (MinGW).

## Building and Running

### Prerequisites
- CMake 3.16+
- C99-compliant compiler (Clang, GCC, or MSVC)
- *(Optional)* FFMpeg libraries ('libavformat', 'libavcodec', 'libavutil', 'libswscale') for direct video decoding.

### Build
```bash
cmake -B build
cmake --build build
```

### Run Unit Tests
```bash
ctest --test-dir build --output-on-failure
# or run directly
./build/unit_tests
```

### CLI Usage
```bash
# 1. Run synthetic test scene (no inputs required)
./build/PALETTE

# 2. Extract a 5-color palette from an image and export swatches
./build/PALETTE -i photo.ppm -k 5 -0 swatch.ppm

# 3. Extract a palette from a video frame
./build/PALETTE -i video.mp4 -k 6

# 4 Use CIELAB color space instead of Oklab
./build/PALETTE -i photo.bmp -k 5 --cielab
```

### Generate Testing Images
(Requires numpy)
```python
python generate_test_images.py
```




