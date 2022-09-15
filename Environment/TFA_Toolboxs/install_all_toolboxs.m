%INSTALL_EMD.M install the EMD package
% add the following locations to Matlab's path:
% <EMD_BASEDIR>
% <EMD_BASEDIR>/emd
% <EMD_BASEDIR>/utils
% <EMD_BASEDIR>/examples/NSIP2003
% <EMD_BASEDIR>/examples/SPL2007
%
% where <EMD_BASEDIR> is the directory containing install_emd.m
%
% and compiles the C codes
%A
% IMPORTANT: After running INSTALL_EMD you must run the "savepath" command to save the installation
% but be careful that if you previously removed parts of the path (using e.g. the "rmpath" command) 
% these will be permanently removed after you run "savepath"

function install_all_toolboxs

base_dir = fileparts(which('install_all_toolboxs'));
addpath(genpath([base_dir,'/EMD']));
addpath(genpath([base_dir,'/TFTB']));
addpath(genpath([base_dir,'/others']));
addpath([base_dir,'/package_emd']);
install_emd
savepath

disp('EMD/TFTB/package_emd工具箱安装完成')
