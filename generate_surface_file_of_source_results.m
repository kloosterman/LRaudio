%% Generate surface file of source results
%
% The codes below can be used in combination with source estimations and / or
% statistical analysis to generate the surface file of the results. 
% The surface file is saved on disk and can be visualized using HCP Connectome Workbench. 
%                                                                         
% Requires:
%         - From BALSA parcellation (https://balsa.wustl.edu/WN56) CIFTI File:
%               Q1-Q6_RelatedParcellation210.CorticalAreas_dil_Colors.32k_fs_LR.dlabel.nii
%         - FieldTrip (fieldtrip-master from Github; ft_read_cifti.m, ft_write_cifti.m, gifti)
%           
% Mohsen Alavash Feb 2022                                               
% 
%% (1) Configuration
%  -----------------
clear
clc
close all

disp('   *** Generating surface file of source results ... ')

HomeDIR         = 'INSERT PATH';
DataDIR         = 'INSERT PATH';
GiftiDIR        = 'INSERT PATH, FOR EXAMPLE: D:/fieldtrip-master/external/gifti';
FieldTripDIR    = 'INSERT PATH, FOR EXAMPLE: D:/fieldtrip-master';
ParcellationDIR = 'INSERT PATH, FOR EXAMPLE: D:/my project/Parcellation';
DestinDIR       = 'INSERT PATH';

addpath(genpath(HomeDIR))
addpath(GiftiDIR)
addpath(FieldTripDIR)
ft_defaults

% Information about the parcellation
no_of_nodes = 362;

% IDs of parcels across right hemisphere (RH)
R_Occi  = [181,182,183,184,185,186,187,193,195,196,197,198,199,200,201,202,203,301,318,320,321,322,323,326,332,333,334,336,337,338,339,340,343]; 
R_Pari  = [188,189,194,207,209,210,211,212,214,213,215,216,217,218,219,220,221,222,223,225,226,227,228,229,230,231,232,233,234,235,236,275,276,280,281,282,295,296,297,324,325,327,328,329,330,331,341,342,348];
R_Temp  = [204,205,208,283,284,285,286,287,290,292,298,299,300,302,303,304,305,306,307,308,309,310,311,312,313,314,315,316,317,319,335,347,352,353,354,355,356,357,358];
R_Front = [190,191,192,206,224,237,238,239,240,241,242,243,244,245,246,247,248,249,250,251,252,253,254,255,256,257,258,259,260,261,262,263,264,265,266,267,268,269,270,271,272,273,274,277,278,279,288,289,291,293,294,344,345,346,349,350,351,359,360];

% LH
L_Occi  = R_Occi -180; 
L_Pari  = L_Pari -180;
L_Temp  = L_Temp -180;
L_Front = L_Front -180;

%% (2) Source analysis and / or statistical analysis
%  -------------------------------------------------
% In section 3 it is assumed that the final results of the source
% analysis implemented in section 2 is saved in the variable S (size(S): no_of_nodes x 1)

%% (3) Surface visualization of source results
% --------------------------------------------
outdir      = DestinDIR;
outfilename = 'RESULTS.32k_fs_LR.dlabel.nii'; % Only adapt RESULTS as you wish, keep the file suffix as it is.

addpath(FieldTripDIR)
ft_defaults

disp('   *** Surface generation of source results ... ')
disp('   *** Reading the CIFTI file of parcellation labels ... ')

cd(ParcellationDIR)
incifti = ft_read_cifti('Q1-Q6_RelatedParcellation210.CorticalAreas_dil_Colors.32k_fs_LR.dlabel.nii');

disp('   *** DONE.')

outcifti = incifti;
outcifti.x1 = nan(length(incifti.x1), 1);

% LH effects
for n = 1:180
    i = find(incifti.x1==n);
    outcifti.x1(i) = S(n);
end

% RH effects 
for n = 181:360
    i = find(incifti.x1==n);
    outcifti.x1(i) = S(n);
end

% Flip LR
fsL = find(incifti.brainstructure==1);
fsR = find(incifti.brainstructure==2);

outcifti.brainstructure(fsL) = 2;
outcifti.brainstructure(fsR) = 1;

disp('   *** Writing the CIFTI file of source results ... ')

cd(outdir)
ft_write_cifti(outfilename, outcifti, 'parameter', 'x1')   

disp('   *** DONE.')

