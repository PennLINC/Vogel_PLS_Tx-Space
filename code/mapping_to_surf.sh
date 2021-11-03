### NOTE!
### In order to use this script, you will need to have connectome workbench installed
### Change the very first command to of each line to your local wb_command or alias.

#C1
/Users/jacobv/Downloads/workbench/bin_macosx64/wb_command -volume-to-surface-mapping ../outputs/c1_xp_image.nii.gz ../surf/S1200.L.midthickness_MSMAll.32k_fs_LR.surf.gii ../surf/PLS_Volume_C1_Smth0mm_gc.func.gii -enclosing

/Users/jacobv/Downloads/workbench/bin_macosx64/wb_command -metric-smoothing -fix-zeros ../surf/S1200.L.midthickness_MSMAll.32k_fs_LR.surf.gii ../surf/PLS_Volume_C1_Smth0mm_gc.func.gii 5 ../surf/PLS_Volume_C1_Smth0mm_surf5mm_gc.func.gii

#C2
/Users/jacobv/Downloads/workbench/bin_macosx64/wb_command -volume-to-surface-mapping ../outputs/c2_xp_image.nii.gz ../surf/S1200.L.midthickness_MSMAll.32k_fs_LR.surf.gii ../surf/PLS_Volume_C2_Smth0mm_gc.func.gii -enclosing

/Users/jacobv/Downloads/workbench/bin_macosx64/wb_command -metric-smoothing -fix-zeros ../surf/S1200.L.midthickness_MSMAll.32k_fs_LR.surf.gii ../surf/PLS_Volume_C2_Smth0mm_gc.func.gii 5 ../surf/PLS_Volume_C2_Smth0mm_surf5mm_gc.func.gii

#C3
/Users/jacobv/Downloads/workbench/bin_macosx64/wb_command -volume-to-surface-mapping ../outputs/c3_xp_image.nii.gz ../surf/S1200.L.midthickness_MSMAll.32k_fs_LR.surf.gii ../surf/PLS_Volume_C3_Smth0mm_gc.func.gii -enclosing

/Users/jacobv/Downloads/workbench/bin_macosx64/wb_command -metric-smoothing -fix-zeros ../surf/S1200.L.midthickness_MSMAll.32k_fs_LR.surf.gii ../surf/PLS_Volume_C3_Smth0mm_gc.func.gii 5 ../surf/PLS_Volume_C3_Smth0mm_surf5mm_gc.func.gii