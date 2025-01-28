# Guide 

This package implements methods for diagnosing and correcting spatial unit roots developed by Müller and Watson (2024). A practical guide to these methods and the Stata implementation can be found in Becker, [Boll and Voth (2025)](https://warwick.ac.uk/fac/soc/economics/research/workingpapers/2025/twerp_1541-_becker.pdf).

When using this code, please cite [Becker, Boll and Voth (2025)](https://warwick.ac.uk/fac/soc/economics/research/workingpapers/2025/twerp_1541-_becker.pdf):

```bash

@TechReport{becker2025,
  author={Becker, Sascha O. and Boll, P. David and Voth, Hans-Joachim},
  title={Spatial Unit Roots in Regressions: A Practitioner's Guide and a Stata Package},
  year=2025,
  institution={University of Warwick, Department of Economics},
  type={The Warwick Economics Research Paper Series (TWERPS)},
  url={https://warwick.ac.uk/fac/soc/economics/research/workingpapers/2025/twerp_1541-_becker.pdf},
  number={1541},
}

```

### Installation

To install the package, use the following command inside your R script:

```bash
install.packages("spur")
library(spur)
```