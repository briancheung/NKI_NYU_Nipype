#!/usr/bin/env bash

##########################################################################################################################
## SCRIPT TO PREPROCESS THE FUNCTIONAL SCAN
## parameters are passed from 0_preprocess.sh
##
## Written by the Underpants Gnomes (a.k.a Clare Kelly, Zarrar Shehzad, Maarten Mennes & Michael Milham)
## for more information see www.nitrc.org/projects/fcon_1000
##
##########################################################################################################################

## subject
subject=0010052
## analysisdirectory/subject
dir=/home/sharad/Resources/0010052
## name of the resting-state scan
rest=rest
## first timepoint (remember timepoint numbering starts from 0)
TRstart=0
## last timepoint
TRend=175
## TR
TR=2.0
## session
session=session_1
##routine
routine=rest_1

## set your desired spatial smoothing FWHM - we use 6 (acquisition voxel size is 3x3x4mm)
FWHM=6.0
sigma=`echo "scale=10 ; ${FWHM}/2.3548" | bc`
echo "spatial smoothing FWHM = $FWHM ( should be 1.5 ~ 2 times of voxel size, e.g. FWHM=6 for 3x3x4mm voxel )"

## Set high pass and low pass cutoffs for temporal filtering
shift; hp=0.005
shift; lp=0.1

## directory setup
func_dir=${dir}/${session}/$routine
mkdir $func_dir
### copy needed file
cp ${dir}/${session}/*.* $func_dir

##########################################################################################################################
##---START OF SCRIPT----------------------------------------------------------------------------------------------------##
##########################################################################################################################

echo ---------------------------------------
echo !!!! PREPROCESSING FUNCTIONAL SCAN: ${subject} ${session} ${routine} !!!!
echo ---------------------------------------

cwd=$( pwd )
cd ${func_dir}

if [ -e "${rest}_pp_mask.nii.gz" ]
then
	echo -----------------------------------------------------------------------------------------
	echo 		!!! ${subject} initial preprocessing already complete, SKIPPING !!!
	echo If you want to rerun functional preprocessing, DELETE ${func_dir}/${rest}_pp.nii.gz FIRST
	echo
else

	## 1. Dropping first # TRS
	echo "Dropping first TRs"; echo "[${TRstart}..${TRend}]"
	3dcalc -a ${rest}.nii.gz[${TRstart}..${TRend}] -expr 'a' -prefix ${rest}_dr.nii.gz

	##2. Deoblique
	echo "Deobliquing ${subject}"
	3drefit -deoblique ${rest}_dr.nii.gz

	##3. Reorient into fsl friendly space (what AFNI calls RPI)
	echo "Reorienting ${subject}"
	3dresample -orient RPI -inset ${rest}_dr.nii.gz -prefix ${rest}_ro.nii.gz

	
		##4. Motion correct to average of timeseries ## Yang - changed Feb 16, 2012
		echo "Motion correcting ${subject}"
		# Initial reg
		3dTstat -mean -prefix ${rest}_ro_mean.nii.gz ${rest}_ro.nii.gz 
		3dvolreg -Fourier -twopass -base ${rest}_ro_mean.nii.gz -zpad 4 \
			-prefix ${rest}_mc_init.nii.gz -1Dfile ${rest}_mc_init.1D \
			-maxdisp1D ${rest}_maxdisp_init.1D ${rest}_ro.nii.gz
		$ 2nd step reg
		3dTstat -mean -prefix ${rest}_mc_init_mean.nii.gz ${rest}_mc_init.nii.gz 
		3dvolreg -Fourier -twopass -base ${rest}_mc_init_mean.nii.gz -zpad 4 \
			-prefix ${rest}_mc.nii.gz -1Dfile ${rest}_mc.1D \
			-maxdisp1D ${rest}_maxdisp.1D ${rest}_ro.nii.gz



	## 5. Remove skull/edge detect
	echo "Skull stripping ${subject}"
	3dAutomask -prefix ${rest}_mask.nii.gz -dilate 1 ${rest}_mc.nii.gz
	3dcalc -a ${rest}_mc.nii.gz -b ${rest}_mask.nii.gz -expr 'a*b' -prefix ${rest}_ss.nii.gz

	##6. Get eighth image for use in registration ## Yang - changed Feb 16, 2012
	echo "Getting example_func for registration for ${subject}"
	#3dcalc -a ${rest}_ss.nii.gz[4] -expr 'a' -prefix example_func.nii.gz
	3dTstat -mean -prefix example_func.nii.gz ${rest}_ss.gz 
		

	##6.5 (ADDED 12th August 2011) Despiking ## Yang - changed Feb 16, 2012, no despiking??
	##echo "Despiking ${subject}"
	##3dDespike -prefix  ${rest}_ds.nii.gz ${rest}_ss.nii.gz
	
	##7. Spatial smoothing  ## -commented out ## Yang - changed Feb 16, 2012
	##echo "Smoothing ${subject}"
	##fslmaths ${rest}_ds.nii.gz -kernel gauss ${sigma} -fmean -mas ${rest}_mask.nii.gz ${rest}_sm.nii.gz

	##8. Grandmean scaling
	echo "Grand-mean scaling ${subject}"
	#fslmaths ${rest}_sm.nii.gz -ing 10000 ${rest}_gms.nii.gz -odt float
	fslmaths ${rest}_ss.nii.gz -ing 10000 ${rest}_gms.nii.gz -odt float

	##11.Create Mask
	cp ${rest}_gms.nii.gz ${rest}_pp.nii.gz ## added by Yang, to keep the naming convention, Feb 16, 2012
	echo "Generating mask of preprocessed data for ${subject}"
	fslmaths ${rest}_pp.nii.gz -Tmin -bin ${rest}_pp_mask.nii.gz -odt char

fi

gzip *.nii

echo ---------------------------------------
echo !!!! PREPROCESSING FUNCTIONAL SCAN Done: ${subject} ${session} ${routine} !!!!
echo ---------------------------------------

pwd

#echo $( pwd ) >> ${cwd}/ppcount.log
ls *_pp* | tee -a ${cwd}/ppcount.log



cd ${cwd}
