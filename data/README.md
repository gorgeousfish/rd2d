# rd2d Package Data Files

## data_rd2d.csv

Synthetic location-based boundary RD dataset (N = 20,000).

| Variable | Type | Description |
|----------|------|-------------|
| x.1 | numeric | First running variable (horizontal coordinate) |
| x.2 | numeric | Second running variable (vertical coordinate) |
| t | logical | Treatment indicator (TRUE if unit is in treated region) |
| y | numeric | Outcome variable |

**Treatment assignment**: t = TRUE when unit lies in the treated region
defined by the boundary geometry.

**Data generating process**: Synthetic data with known treatment effect
for package testing and demonstration.

## D.csv

Signed-distance matrix (40 × 40) for the distance-based RD commands.
Each row corresponds to one evaluation point in `eval.csv`.
Each column represents a signed distance variable for a boundary segment.

## eval.csv

Evaluation point coordinates (40 × 2) for boundary RD estimation.

| Variable | Type | Description |
|----------|------|-------------|
| x.1 | numeric | First coordinate of boundary evaluation point |
| x.2 | numeric | Second coordinate of boundary evaluation point |

These 40 points define locations along the treatment boundary where
local polynomial RD estimates are computed.

## Usage

```stata
* Location-based RD
import delimited "data_rd2d.csv", clear
rd2d y x1 x2 t, at(0 25) p(1) kernel(triangular)

* Distance-based RD
import delimited "D.csv", clear
rd2d_dist y dist*, p(1) kernel(triangular)
```

## Provenance

All data files are synthetically generated for package testing
and demonstration. No real-world data is included.
