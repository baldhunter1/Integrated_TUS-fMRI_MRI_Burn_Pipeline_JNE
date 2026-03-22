# MRI Burn Pipeline

Script for the MRI burn pipeline for the publication "Integrated Neuronavigation and Acoustic-Thermal Modeling for Simultaneous Transcranial Ultrasound Stimulation and Functional MRI."

## Description

This script performs rigid-body registration of a high-resolution T1w image and a target binary mask into a rapid/scout T1w space using FSL FLIRT, then burns the registered mask into DICOM images as a white intensity overlay (dot). The script runs in a graphical user interface (GUI).

## Requirements

- Python 3
- [FSL](https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/) (FLIRT must be available in PATH)
- Python packages: `nibabel`, `pydicom`, `numpy`, `tkinter`
- [dcm2niix](https://github.com/rordenlab/dcm2niix)

## Usage

Run the script:

```bash
python mri_burn_pipeline.py
```

The GUI will prompt for the following inputs:

| Input | Description |
|---|---|
| Standard T1 (NIfTI) | High-resolution/native T1w image |
| Rapid/Scout T1 (NIfTI) | Low-resolution T1w image used as registration reference |
| Native space Mask | Binary mask of the target ROI in native T1w space |
| Output folder | Directory where registration outputs will be saved |
| DICOM IN folder | Folder containing the original DICOM files |
| DICOM OUT folder | Folder where burned DICOM files will be saved |

All images must be in NIfTI format (`.nii` or `.nii.gz`).

## Outputs

- `registered_t1w.nii.gz` — T1w image registered to rapid/scout space
- `registered_t1w_mask.nii.gz` — Binary mask registered to rapid/scout space
- `std2rapid.mat` — FLIRT transformation matrix
- Burned DICOM files saved to the specified DICOM OUT folder

## Acknowledgements

This script is adapted from the MRI burn pipeline code provided by Dr. Samuel Pichardo from the University of Calgary.
