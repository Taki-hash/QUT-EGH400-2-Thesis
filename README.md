# Buck Converter THD Analysis for EIS Applications

This repository supports the thesis titled **"Linearity and General Analysis of Sinusoidal Injection Locations in Buck Converters for Electrochemical Impedance Spectroscopy"**, completed as part of EGH400-2 at Queensland University of Technology (QUT) located in Brisbane, QLD.

## 🧪 Project Summary

This project analyzes how sinusoidal signal injections at different nodes in a buck converter affect Total Harmonic Distortion (THD) and frequency response characteristics. The goal is to determine the most suitable injection location for EIS measurements on DC-DC coupled PEM electrolysis systems.

Key contributions:
- Implementation of a nonlinear buck converter model.
- THD calculation and summation across injection signal with respect to two relevant fundemental frequencies.
- Comparison of input, output, and PWM injection locations.
- Analysis of the trade-off between bandwidth and distortion.

## 🔧 Requirements

- Used MATLAB Version R2024b
- Control System Toolbox (for Bode analysis)
- Signal Processing Toolbox (for FFT)

## 📈 How to Use

1. Clone or download the repository.
2. Run the desired analysis script (`NonlinearAnalysis_IN.m`, etc.).
3. Use the plots and THD results to compare injection points.
4. All scripts are commented and modular for modification or integration.

## 📜 Licensing

This work is part of an academic research thesis and is provided for educational use only.

## 📚 Citation

Evans, K. T. (2025). *Investigation of Harmonic Distortion and Frequency Response in Buck Converters for Electrochemical Impedance Spectroscopy (EIS)*. Undergraduate Thesis, Queensland University of Technology.

@thesis{evans2025buck,
  author       = {Kyle T. Evans},
  title        = {Investigation of Harmonic Distortion and Frequency Response in Buck Converters for Electrochemical Impedance Spectroscopy (EIS)},
  year         = {2025},
  school       = {Queensland University of Technology},
  type         = {Undergraduate Thesis},
  url          = {https://github.com/[your-username]/[your-repo-name]}
}


