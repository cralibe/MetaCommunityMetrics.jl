<table>
  <tr>
    <td><img src="docs/src/assets/logo.png" alt="Logo" width="100" height="100"></td>
    <td><h1>MetaCommunityMetrics</h1></td>
  </tr>
</table>



`MetaCommunityMetrics` is a collection of tools and utilities for analyzing meta-communities in Julia. The current version is compatible with `julia` version 1.9.3.

## An Overview
This package is a comprehensive toolkit designed to characterize the spatiotemporal structure and dynamics of a metacommunity—a network of communities linked by the dispersal of multiple, interacting species, each with unique niche breadths. It includes functions to calculate a range of specific metrics, which have been previously implemented in R and proven valuable for metacommunity analysis. 

However, they come with high computational costs, especially for large species community datasets. To address this issue, MetaCommunityMetrics.jl was developed in Julia, a programming language known for its efficiency in handling computationally intensive tasks. This implementation significantly improves the efficiency of calculating these metrics, making it a powerful tool for metacommunity analysis. 

These metrics include:
- Averaged beta diversity decomposition in space/time: total diversity, species replacement (turnover), and richness differences for both presence-absence and abundance data
- Dispersal-niche continuum index to evaluate the degree to which communities are influenced by dispersal processes and niche breadth
- Niche overlap indices to determine the extent of niche sharing among species within the metacommunity
- The proportion of habitat patches occupied by each species
- The variability of community composition across different spatial and temporal scales
- Niche hypervolume measurements (individual species, average, and between-species dissimilarities)


## Getting Started

### Installation

To install MetaCommunityMetrics, use the following command:

```julia
using Pkg

Pkg.add("MetaCommunityMetrics")

using MetaCommunityMetrics
```
### Sample Data
To assess the sample data, use the following command:
```julia
using MetaCommunityMetrics

load_sample_data()
```

### Example
```@jildoctest
julia> using MetaCommunityMetrics

julia> load_sample_data()
48735×12 DataFrame
   Row │ Year   Month  Day    Sampling_date_order  plot   Species  Abundance  Presence  Latitude  Longitude  normalized_temperature  normalized_precipitation 
       │ Int64  Int64  Int64  Int64                Int64  String3  Int64      Int64     Float64   Float64    Float64                 Float64                  
───────┼──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
     1 │  2010      1     16                    1      1  BA               0         0      35.0     -110.0               0.80987                  -0.290381
     2 │  2010      1     16                    1      2  BA               0         0      35.0     -109.5              -1.12523                   0.750317
     3 │  2010      1     16                    1      8  BA               0         0      35.5     -109.5              -1.10775                  -1.87583
     4 │  2010      1     16                    1      9  BA               0         0      35.5     -109.0              -0.418417                  0.0964911
     5 │  2010      1     16                    1     11  BA               0         0      35.5     -108.0               0.287892                 -0.0272079
   ⋮   │   ⋮      ⋮      ⋮             ⋮             ⋮       ⋮         ⋮         ⋮         ⋮          ⋮                ⋮                        ⋮
 48731 │  2023      3     21                  117      9  SH               0         0      35.5     -109.0               0.143276                  2.37981
 48732 │  2023      3     21                  117     10  SH               0         0      35.5     -108.5               0.148338                  1.4683
 48733 │  2023      3     21                  117     12  SH               1         1      35.5     -107.5               1.01169                  -0.485298
 48734 │  2023      3     21                  117     16  SH               0         0      36.0     -108.5               0.0284359                -0.392446
 48735 │  2023      3     21                  117     23  SH               0         0      36.5     -108.0               0.170814                  1.14892
                                                                                                                                           48725 rows omitted
```

## Acknowledgment
This package includes translations and adaptations of functions from the R packages `adespatial` (licensed under GPL-3), `vegan` (licensed under GPL-2 or later), and `DNCImper` (licensed under GPL-3). The original packages and their documentation are available at:

- `adespatial`: [https://cran.r-project.org/web/packages/adespatial/index.html](https://cran.r-project.org/web/packages/adespatial/index.html)
- `vegan`: [https://cran.r-project.org/web/packages/vegan/index.html](https://cran.r-project.org/web/packages/vegan/index.html)
- `DNCImper`: [https://github.com/Corentin-Gibert-Paleontology/DNCImper](https://github.com/Corentin-Gibert-Paleontology/DNCImper)


Please refer to these sources for full details on the original implementations and licenses.


[![Build Status](https://github.com/cralibe/MetaCommunityMetrics.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/cralibe/MetaCommunityMetrics.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://cralibe.github.io/MetaCommunityMetrics.jl/)
[![codecov](https://codecov.io/github/cralibe/MetaCommunityMetrics.jl/graph/badge.svg?token=OKUWBS8R7U)](https://codecov.io/github/cralibe/MetaCommunityMetrics.jl)

## License
This project is licensed under the terms of the GNU General Public License v3.0. See the LICENSE file for more details.
