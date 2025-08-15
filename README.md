# MetFamily <a href="https://ipb-halle.github.io/MetFamily/index.html"><img src="inst/MetFamily/www/img/Metfamily.gif" align="right" height="120" /></a>

[![build](https://github.com/ipb-halle/MetFamily/actions/workflows/build.yml/badge.svg)](https://github.com/ipb-halle/MetFamily/actions/workflows/build.yml)

## Overview

Understanding metabolism is fundamental in biomedical and plant research and the
identification and quantification of thousands of metabolites by mass
spectrometry in modern metabolomics is a prerequisite for elucidating this area.
However, the identification of metabolites is a major bottleneck in traditional
approaches hampering advances. Here, we present a novel approach for the
untargeted discovery of metabolite families offering a bird's eye view of
metabolic regulation in comparative metabolomics. We implemented the presented
methodology in the easy-to-use web application MetFamily to enable the analysis
of comprehensive metabolomics studies for all researchers worldwide.

## Using MetFamily

The MetFamily web application is available at
https://msbi.ipb-halle.de/MetFamilyDevel/.

It runs inside a Kubernetes cluster at the 
[Leibniz Institute of Plant Biochemistry (IPB)](https://www.ipb-halle.de/en/).

The release version described in [Treutler et al., 2016](https://pubs.acs.org/doi/10.1021/acs.analchem.6b01569) is hosted at https://msbi.ipb-halle.de/MetFamily/.


### Local installation of the R Package

The MetFamily R package can be installed using:
```
remotes::install_github("ipb-halle/MetFamily")
```

After installation, you can run the interactive application locally:
```
library(MetFamily)
runMetFamily()
```

Or use an R script to process your analysis. See our vignette [Discovering regulated Metabolite Families](https://ipb-halle.github.io/MetFamily/articles/discoveringregulatedmetabolitefamilies.html).


## Development

MetFamily is under active development. Let us know if you run into issues.

Documentation for developers is stored in the `dev` folder.

