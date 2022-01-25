#!/bin/bash
# Applies downsampling, masking, mean division normalization, and smoothing
slideDir=$1
# downsample percentage in total pixels
dp=$2
# mask threshold
maskThr=$3
# smoothing kernel
sigma=$4

# downsample all images in directory
echo "downsampling"
for img in $(ls $slideDir/*region*[0-9]{.jpg,.png} 2>/dev/null); do
  oimg=$(dirname $img)/$(basename $img | sed "s+\.[a-z][a-z][a-z]+_downsampled.tif+g")
  c2d $img -resample $dp% -o $oimg
  #c2d $oimg -info-full
done

# loop through regions
regionSuffix=$(ls $slideDir/DAPI_ADJ_region*.jpg | sed "s+.*_region_+region_+g" | sed "s+\.[a-z][a-z][a-z]++g")
for region in $regionSuffix; do
  # create mask
  echo "creating mask for $region"
  dapi=$slideDir/DAPI_ADJ_${region}_downsampled.tif
  tissueMask=$slideDir/mask_${region}_downsampled.tif
  c2d $dapi -threshold $maskThr inf 1 0 -o $tissueMask

  # divide each image by the mean within the DAPI mask
  for img in $(ls $slideDir/*${region}_downsampled.tif);do
    echo "  normalizing and smoothing $(basename $img)"
    mean=$(c2d $img $tissueMask -lstat | awk 'FNR == 3 {print $2}')
    meaninv=$(echo "1/$mean" | bc -l)
    oimg=$(echo $img | sed "s+\.tif+_normed.tif+g")
    c2d $img -scale $meaninv -o $oimg

    # Apply smoothing
    # -ad is anisotropic diffusion
    ooimg=$(echo $oimg | sed "s+\.tif+_smoothed.tif+g")
    c2d $oimg -smooth "${sigma}x${sigma}vox" -o $ooimg
  done

done
