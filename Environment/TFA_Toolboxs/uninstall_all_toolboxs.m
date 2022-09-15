% UNINSTALL_EMD.M uninstall the EMD package
% reverts the modifications made by INSTALL_EMD:
% removes the following locations from Matlab's path:
% <EMD_BASEDIR>
% <EMD_BASEDIR>/utils
% <EMD_BASEDIR>/examples/NSIP2003
% <EMD_BASEDIR>/examples/SPL2007
%
% where <EMD_BASEDIR> is the directory containing uninstall_emd.m
%
% UNINSTALL_EMD also interactively proposes to remove the files from the hard drive
%
% IMPORTANT: After running UNINSTALL_EMD you must run the "savepath" command to save the modifications made 
% to Matlab's path but be careful that if you previously removed parts of the path (using e.g. the "rmpath" command) 
% these will be permanently removed after you run "savepath"
function uninstall_all_toolboxs
base_dir = fileparts(which('uninstall_all_toolboxs'));
rmpath(base_dir)
rmpath(genpath([base_dir,'/EMD']))
rmpath(genpath([base_dir,'/TFTB']))
rmpath(genpath([base_dir,'/others']))
%uninstall_emd
rmpath(genpath([base_dir,'/package_emd']))

disp('EMD/TFTB/package_emd工具箱卸载完成')
