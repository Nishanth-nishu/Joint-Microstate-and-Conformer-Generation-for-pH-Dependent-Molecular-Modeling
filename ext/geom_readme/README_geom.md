This README file was generated on 2022-02-10 by SIMON AXELROD.


GENERAL INFORMATION

1. Title of Dataset: GEOM

2. Author Information


	A. First author
	
		Name: Simon Axelrod
		Institution: Harvard University
		Address: 12 Oxford Street, Cambridge, MA, USA
		Email: simonaxelrod83@gmail.com

	B. Second author
		
		Name: Rafael Gomez-Bombarelli
		Institution: Massachusetts Institute of Technology
		Address: 77 Massachusetts Avenue, Cambridge, MA, USA
		Email: rafagb@mit.edu


3. Date of data generation: May 2020

4. Information about funding sources that supported the collection of the data: XSEDE COVID-19 HPC Consortium, project CHE200039. Computational resources provided by NASA Advanced Supercomputing (NAS) Division and LBNL National Energy Research Scientific Computing Center (NERSC), MIT Engaging
cluster, Harvard Cannon cluster, and MIT Lincoln Lab Supercloud.



SHARING/ACCESS INFORMATION

1. Licenses/restrictions placed on the data: CC0 1.0 Universal.

2. Links to publications that cite or use the data: https://arxiv.org/abs/2006.05531, https://arxiv.org/abs/2012.08452

3. Links to other publicly accessible locations of the data: https://www.dropbox.com/sh/1aptf9fi8kyrzg6/AABQ4F7dpl4tQ_pGCf2izd7Ca?dl=0

4. Links/relationships to ancillary data sets: None

5. Was data derived from another source? No

6. Recommended citation for this dataset: Bibtex: @article{axelrod2020geom,
  title={GEOM: Energy-annotated molecular conformations for property prediction and molecular generation},
  author={Axelrod, Simon and Gomez-Bombarelli, Rafael},
  journal={arXiv preprint arXiv:2006.05531},
  year={2020}
}


DATA & FILE OVERVIEW

1. File list
	- There are two tarred messagepack files with just geometries:
		- `drugs_crude.msgpack.tar.gz`
		- `qm9_crude.msgpack.tar.gz`
	- There are two tarred messagepack files with additional graph-based features:
		- `drugs_featurized.msgpack.tar.gz`
		- `qm9_featurized.msgpack.tar.gz`
	- There are six tarred folder of pickle files:
		- `rdkit_folder.tar.gz`
		- `bace_water.tar.gz`
		- `censo.tar.gz`
		- `molecule_net.tar.gz`
		- `censo_hess_and_orbs.tar.gz`
		- `crest_hess_and_orbs.tar.gz`
		
	- The tarred files can be decompressed by running `tar -xzvf <tarred_file_name>`. Decompressing the `.msgpack.tar.gz` files will give you messagepack files, while decompressing `rdkit_folder.tar.gz` will give you a folder.
	- The first four pickle files represent the geometries as RDKit mol objects. Each species has its own pickle file with all of its
	conformers. The last two pickle files contain extra data on the computed Hessians and orbitals for CENSO and CREST geometries in the BACE dataset.
	- The messagepack files are language-agnostic, whereas the pickle files use Python, but the pickle files are much easier to handle. This is because RDKit can perform all sorts
	of analysis on Mol objects without the user having to write any of it. 
	- **For this reason we will not be updating any of the messagepack files with new data; only the RDKit folder will be updated. If you want to get the molecules in non-Python form, please contact us and we will be happy to help out.**


2. Data updates

	i) New drug-like molecules (Feb. 1, 2021). We have updated the GEOM dataset since our paper was first posted on the ArXiv, adding about 13,000 new drug-like molecules, including about 6,000 with SARS-CoV-2 data. 
	
	ii) MoleculeNet (Feb. 9, 2022). We have added over 16,000 species from the MoleculeNet dataset. These species have data related to physical chemistry, physiology, and biophysics. These results can be found with the rest of the data here. The data contains CREST ensembles for the species, together with Hessian data and high-accuracy DFT results for the BACE subset. 
	
	- Checksums are used to confirm that your version of the data matches our version, by computing the message digest of a file.
	- MD5 checksums are automatically generated for all files on Dataverse, and they are shown two lines below the file name. You can get the MD5 message digest of your file and compare to the value shown on Dataverse to make sure you have the right version of the file.
	- To get the MD5 message digest of a file, run `md5sum <your_file>` on linux, `certutil -hashfile <your_file> md5` on Windows, or `md5 <your_file>` on mac. This may take a few minutes as the files are quite big.
	
	
3. Loading the data
	- For tutorials on loading the data, please use our [GEOM github repository](https://github.com/learningmatter-mit/geom)
	

