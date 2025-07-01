# Integrating multiplexing into confineable gene drives effectively overrides resistance in Anopheles stephensi

**Authors:**  
Mireia Larrosa Godall¹³⁴, Lewis Shackleford¹²³, Matthew P. Edgington¹²³, Philip T. Leftwich⁵, James C. Y. Luk¹², Joshua Southworth³, Stewart Rosell³, Jake Creasey¹², Jack Aked¹², Katherine Nevard³, Alexander Dodds¹², Morgan Mckee¹², Eunice Adedeji¹², Estela Gonzalez³⁶, Joshua X.D. Ang¹²³, Michelle A. E. Anderson¹²³, Luke Alphey¹²³

- ¹ Department of Biology, University of York, UK  
- ² York Biomedical Research Institute, University of York, UK  
- ³ Arthropod Genetics, The Pirbright Institute, UK  
- ⁴ The Francis Crick Institute, UK  
- ⁵ School of Biological Sciences, University of East Anglia, UK  
- ⁶ Animal and Plant Health Agency, UK  

**Keywords:** mosquito, CRISPR-Cas9, malaria

---

## Summary

This repository contains all data, R scripts, and supplementary materials for our study on multiplexed CRISPR/Cas9 gene drives in Anopheles stephensi. Our research demonstrates that multiplexing with multiple sgRNAs can overcome resistance at the target site, a major challenge in gene drive deployment. We provide both experimental data and modeling results to support our findings. All analyses are reproducible using the included code and environment management via `renv`.

---

## Repository Structure

```text
.
├── data/             # Raw and processed data files
├── scripts/          # R scripts for analysis and modeling and CRISPRESSO bash script
├── figures/          # Output figures and plots
├── renv/             # renv environment folder
├── renv.lock         # Lockfile for R package versions
└── README.md         # This file
```

---

## Getting Started

### Prerequisites

- [R](https://www.r-project.org/) (version 4.0 or higher recommended)
- [renv](https://rstudio.github.io/renv/) (for R package and environment management)

### Installation

1. **Clone the repository:**
   ```sh
   git clone https://github.com/Philip-Leftwich/stephensii_multiplexing.git
   cd stephensii_multiplexing
   ```

2. **Restore the project environment using renv:**
   Open R or RStudio in the project directory and run:
   ```r
   renv::restore()
   ```
   This will install all required packages as specified in `renv.lock`.

---

## Usage

1. **Data:**  
   - All raw and processed datasets are in the `data/` directory.

2. **Scripts:**  
   - Open the relevant R scripts in the `scripts/` directory.
   - Each script includes comments and instructions for reproducing analyses and generating figures.
   - Automated CRISPResso2 Pipeline Install [CRISPResso2](https://crispresso.pinellolab.partners.org/) and its dependencies.

3. **Figures:**  
   - Output figures are saved in the `figures/` directory after running scripts.

---

## License

This repository is licensed under the [MIT License](./LICENSE).

---

## Contact

For questions, please contact:  
- michelle.anderson@york.ac.uk  
- luke.alphey@york.ac.uk
- p.leftwich@uea.ac.uk 

---

## Acknowledgements

This project was supported by the University of York, The Pirbright Institute, and collaborators as detailed above.
