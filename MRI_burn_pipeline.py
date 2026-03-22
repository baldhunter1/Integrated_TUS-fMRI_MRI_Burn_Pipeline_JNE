#!/usr/bin/env python3

'''
Script is adapted from the MRI burn pipeline code provided by 
Dr. Samuel Pichardo from the University of Calgary.

This script runs in a GUI and asks the user for a high-resolution/native T1w 
image, a rapid/scout T1w image, a standard-space binary mask (ROI/target), 
an output folder for registration results, and DICOM input/output folders.

It performs a rigid-body registration (FSL FLIRT) of the high-resolution T1w 
and mask into the rapid/scout T1w space, then burns the registered mask into 
the DICOM images as a white intensity overlay.

All images should be in NIfTI format.

Author: Moon Jeong

'''

import subprocess
import sys
import shutil
from pathlib import Path
import tkinter as tk
from tkinter import filedialog, messagebox
import numpy as np
import nibabel
import pydicom
import os
from glob import glob

# FLIRT registration
def run_flirt(std_t1: Path, rapid_t1: Path, mask: Path, outdir: Path):
    outdir.mkdir(parents=True, exist_ok=True)
    mat      = outdir / "std2rapid.mat"
    reg_t1   = outdir / "registered_t1w.nii.gz"
    reg_mask = outdir / "registered_t1w_mask.nii.gz"

    subprocess.check_call([
        "flirt", "-in", std_t1, "-ref", rapid_t1,
        "-omat", mat, "-out", reg_t1,
        "-dof", "6", "-interp", "sinc"
    ])

    subprocess.check_call([
        "flirt", "-in", mask, "-ref", rapid_t1,
        "-applyxfm", "-init", mat,
        "-interp", "nearestneighbour",
        "-out", reg_mask
    ])

    return reg_t1, reg_mask, mat

# showimages class
class showimages(object):
    def __init__(self,NiT1fname, NiftiTractfname, DCMpath,
                 RelativeOverIntensity=0.75,
                 BoundaryTract=0.05,
                 BurnStyle='WhiteIntensity'):

        dcmnix='/opt/homebrew/bin/dcm2niix'
        LDICOM = glob(DCMpath+os.sep+'Z*')
        if len(LDICOM)==0:
            LDICOM = glob(DCMpath+os.sep+'*.dcm')
        lDic=[]
        for f in LDICOM:
            lDic.append(pydicom.dcmread(f))
        lDic.sort(key=lambda x: tuple(x.ImagePositionPatient),reverse=True)
        ret=os.system(dcmnix + ' -w 1 "' + DCMpath + '"' +os.sep)
        if ret!=0:
            raise ValueError('dcm2niix command failed')
        lNifti = glob(DCMpath+os.sep+'*.nii')
        print(lNifti)
        if (len(lNifti)!=1):
            raise ValueError('one nifti file should be at ' + DCMpath+os.sep)
        
        convNifti=nibabel.load(lNifti[0])
        
        self.Nifti=nibabel.load(NiT1fname)
        self.NiftiTract=nibabel.load(NiftiTractfname)
        self.DCM=lDic
        self.RelativeOverIntensity=RelativeOverIntensity
        self.BoundaryTract=BoundaryTract
        print(self.Nifti.affine,'\n',convNifti.affine)
        assert(np.allclose(self.Nifti.affine,self.NiftiTract.affine))
        assert(np.allclose(self.Nifti.affine,convNifti.affine))
        
        self.NiftiData=self.Nifti.get_fdata()
        self.NiftiTractData=self.NiftiTract.get_fdata()
        
        assert(np.all(np.array(self.NiftiData.shape)==np.array(self.NiftiTractData.shape)))
        assert(np.all(np.array(self.NiftiData.shape)==np.array(convNifti.get_fdata().shape)))
        MaxSignalDCM=0
        for im in self.DCM:
            if MaxSignalDCM<im.pixel_array.max():
                MaxSignalDCM=im.pixel_array.max()
        self.MaxSignalDCM=MaxSignalDCM
        self.BurnStyle=BurnStyle
        if BurnStyle=='WhiteIntensity':
            self.GenerateBurnedDatasetWhiteIntensity()
        else:
            self.GenerateBurnedDatasetContour()
        
    def ReorientNifti2DCM(self,implane):
        return np.rot90(implane.T,2)

    def GenerateBurnedNifti(self,destname='out.nii.gz'):
        NewNifti=nibabel.Nifti1Image(self._BData,affine=self.Nifti.affine)
        NewNifti.to_filename(destname)
        
    def GenerateBurnedDatasetWhiteIntensity(self):
        MaxTracto=self.NiftiTractData.max()
        MaxNiftiT1=self.NiftiData.max()
        NTdata=self.NiftiTractData.copy()
        NTdata/=MaxTracto
        NTdata[NTdata<self.BoundaryTract]=0
        self._BData=self.NiftiData+NTdata*MaxNiftiT1*self.RelativeOverIntensity
    
    def GenerateBurnedFromDCM(self,destdir):
        if os.path.isdir(destdir):
            ldcm=glob(destdir+os.sep+'Z*')
            for l in ldcm:
                os.remove(l)
        else:
            os.makedirs(destdir)
        SeriesDescription=''
        MaxTracto=self.NiftiTractData.max()
        for n,l in enumerate(self.DCM):
            BData=l.pixel_array+\
                (self.ReorientNifti2DCM(self.NiftiTractData[n,:,:])/MaxTracto*0.5*self.MaxSignalDCM).astype(np.int16)
            newDicom=l.copy()
            newDicom.PixelData = BData.tobytes()
            
            rootUID='.'.join(l.SOPInstanceUID.split('.')[:-1])+'.'
            newDicom.SOPInstanceUID=pydicom.uid.generate_uid(prefix=rootUID)
            if SeriesDescription=='':
                SeriesDescription=newDicom.SeriesDescription +'-BurnedMoonDOT'
                rootUID='.'.join(l.SeriesInstanceUID.split('.')[:-1])+'.'
                SeriesUID=pydicom.uid.generate_uid(prefix=rootUID)
            
            newDicom.SeriesDescription=SeriesDescription
            newDicom.SeriesInstanceUID=SeriesUID
            newDicom.save_as(destdir+os.sep+'Z%02i' %(n+1))

# GUI App
class App(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("Rigid T1 + Mask Registration + Burn DICOM")
        self.resizable(False, False)

        self.std_t1   = tk.StringVar()
        self.rapid_t1 = tk.StringVar()
        self.mask     = tk.StringVar()
        self.outdir   = tk.StringVar(value=str(Path.cwd() / "reg_output"))
        self.dcm_in   = tk.StringVar(value='/path/to/DICOM/IN/A')
        self.dcm_out  = tk.StringVar(value='/path/to/DICOM/OUT/A')

        fields = [
            ("Standard T1 (NIfTI):",    self.std_t1,   "file"),
            ("Rapid/Scout T1 (NIfTI):", self.rapid_t1, "file"),
            ("Standard-space Mask:",    self.mask,     "file"),
            ("Output folder:",          self.outdir,   "dir"),
            ("DICOM IN folder:",        self.dcm_in,   "dir"),
            ("DICOM OUT folder:",       self.dcm_out,  "dir"),
        ]

        for row, (label, var, kind) in enumerate(fields):
            tk.Label(self, text=label).grid(row=row, column=0, sticky="e", padx=5, pady=2)
            tk.Entry(self, textvariable=var, width=60).grid(row=row, column=1, padx=5, pady=2)
            tk.Button(self, text="…", width=3,
                      command=lambda v=var, k=kind: self._browse(v, k)
                      ).grid(row=row, column=2, padx=5)

        tk.Button(self, text="Run Registration + Burn DICOM",
                  command=self._run).grid(row=len(fields), column=0,
                                          columnspan=3, pady=10, ipadx=40)

    def _browse(self, var, kind):
        path = (filedialog.askopenfilename() if kind == "file"
                else filedialog.askdirectory())
        if path:
            var.set(path)

    def _run(self):
        paths = [self.std_t1, self.rapid_t1, self.mask, self.dcm_in]
        if not all(Path(p.get()).exists() for p in paths):
            messagebox.showerror("Missing input", "Choose all required files before running.")
            return

        try:
            reg_t1, reg_mask, mat = run_flirt(
                Path(self.std_t1.get()),
                Path(self.rapid_t1.get()),
                Path(self.mask.get()),
                Path(self.outdir.get())
            )

            messagebox.showinfo("Registration done", "Now running burn pipeline.")

            A = showimages(
                DCMpath=self.dcm_in.get(),
                NiT1fname=reg_t1,
                NiftiTractfname=reg_mask,
                BurnStyle='WhiteIntensity',
                RelativeOverIntensity=1.0,
                BoundaryTract=0.25
            )

            A.GenerateBurnedFromDCM(self.dcm_out.get())
            messagebox.showinfo("Done", "Burned DICOMs saved!")

        except subprocess.CalledProcessError as e:
            messagebox.showerror("FLIRT error", str(e))
        except Exception as e:
            messagebox.showerror("Error", str(e))

# main 
if __name__ == "__main__":
    if not shutil.which("flirt"):
        sys.exit("Error: FSL FLIRT not found in PATH.")
    App().mainloop()
