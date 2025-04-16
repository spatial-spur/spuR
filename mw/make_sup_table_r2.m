% Read in Results .. 
X_05_1 = NaN(1,n_design);
X_50_1 = NaN(1,n_design);
X_95_1 = NaN(1,n_design);
X_05_5 = NaN(1,n_design);
X_50_5 = NaN(1,n_design);
X_95_5 = NaN(1,n_design);
str1 = ['_k1_' estr];
str5 = ['_k5_' estr];
for i = 1:7
    dstr = char(Design_lst(i));
    fname = [datdir dstr str1 '_5']; 
    tmp = readmatrix(fname);   
    X_05_1(1,i) = tmp;
    fname = [datdir dstr str1 '_50']; 
    tmp = readmatrix(fname);  
    X_50_1(1,i) = tmp;
    fname = [datdir dstr str1 '_95']; 
    tmp = readmatrix(fname); 
    X_95_1(1,i) = tmp;
   
    fname = [datdir dstr str5 '_5']; 
    tmp = readmatrix(fname);
    X_05_5(1,i) = tmp;
    fname = [datdir dstr str5 '_50']; 
    tmp = readmatrix(fname);  
    X_50_5(1,i) = tmp;
    fname = [datdir dstr str5 '_95']; 
    tmp = readmatrix(fname);  
    X_95_5(1,i) = tmp;
end
% Produce Table
outfile_name = [outdir 'um_mc_tab_' file_str '_supplementarymaterial.txt'];
fileID = fopen(outfile_name,'w');

fprintf(fileID,'R2: k = 1 \n');
fprintf(fileID,'DGP; , \n');
for irow = 1:n_design
    nm = char(Design_name(irow));
    fprintf(fileID,[nm ';']);
    tmp = [X_50_1(1,irow) X_05_1(1,irow) X_95_1(1,irow)];
    fprintf(fileID,'%4.3f (%4.3f,%4.3f) \n',tmp);
end

fprintf(fileID,'\n');
fprintf(fileID,'R2: k = 5 \n');
fprintf(fileID,'DGP; , \n');
for irow = 1:n_design
    nm = char(Design_name(irow));
    fprintf(fileID,[nm ';']);
    tmp = [X_50_5(1,irow) X_05_5(1,irow) X_95_5(1,irow)];
    fprintf(fileID,'%4.3f (%4.3f,%4.3f) \n',tmp);
end



      