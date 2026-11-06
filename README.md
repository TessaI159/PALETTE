# PALETTE
**PALETTE**: *Perceived Average of Light in Each Temporal Tile Extracted*

PALETTE extracts representative color palettes from visual data using k-means clustering in perceptually uniform color spaces.

## Current Status ##
- **SoA Migration Complete **: 'struct Color' now stores contiguous channel arrays ('alpha', 'beta', 'gamma') with 32-byte alignment for direct vector SIMD loads.
- **Color Space Conversions**: Bidirectional conversions between linear sRGB, Oklab, and CIELAB with SSE and AVX2 routies where available.
- **Distance Metricx**: Euclidean distance, DIE94, and CIEDE2000 implementations.
- **KMeans Clustering**: Full k-means++ centroid initialization and clustering pipeline on SoA color batches.
- *Cross-Platform**: Build natively on macOS (Apple Silicon / Intel), Linux, and Windows (MinGW).
- **Test Suite**: 100% passing Unity unit tests covering SoA memory lifecycle, conversion roundtrips, distance metrics, averaging, and k-means clustering.

## Building and Running

### Prerequisites
- CMake 3.16+
- C99-compliant compiler (Clang or GCC)
- *(Optional)* FFMpeg libraries ('libavformat', 'libavcodec', 'libavutil', 'libswscale') for video decoding.

### Build and Run Main Demo
'''bash
cmake -B build
cmake --build build
./build/PALETTE
'''


$$$ Run Unit Tests
'''bash
ctest --test-dir build --output-on-failure
# or run directly
./build/unit_tests
'''
