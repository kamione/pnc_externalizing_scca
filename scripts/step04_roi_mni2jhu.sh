#!/bin/bash

ALEseeds='/projects/students/TYWong/2020_pnc_fc_extbeh/masks/ALEseeds'
ALEseedsJHU='/projects/students/TYWong/2020_pnc_fc_extbeh/outputs/ALEseedsJHU'
atlas='/projects/atlas/JHU/MNI/T1/JHU_MNI_SS_T1_ss_reorient_LR.img'
vtk='/projects/atlas/JHU/MNI/T1/PickAtlasLabel/BMAP_mni2JHULR.vtk'


if [ -d "$ALEseedsJHU" ]; then
    rm -rf $ALEseedsJHU;
fi

cp -R $ALEseeds $ALEseedsJHU

cd $ALEseedsJHU

for filename in *.nii
do
    echo $filename
    name="${filename%.nii}"
    echo $name
    
    fslchfiletype NIFTI_PAIR $filename 
    VolumeReorient -i ${name}.img -o ${name}_reorient.img -op YZ #XZU XUY the same
    ApplyLDDMMImg2Img -i ${name}_reorient.img -r $atlas -o ${name}_MNI2JHU.img -deform $vtk -int nearest
 
done
