% Read in Results .. 
X_05_1 = NaN(2,ncol,n_design);
X_50_1 = NaN(2,ncol,n_design);
X_95_1 = NaN(2,ncol,n_design);
X_05_5 = NaN(2,ncol,n_design);
X_50_5 = NaN(2,ncol,n_design);
X_95_5 = NaN(2,ncol,n_design);
str1 = ['_k1_' estr];
str5 = ['_k5_' estr];
for i = 1:n_design
    dstr = char(Design_lst(i));
    fname = [datdir dstr str1 '_5']; 
    tmp = readmatrix(fname);  
    tmp = packr(tmp')';
    X_05_1(:,:,i)=tmp(n_row:n_row+1,1:ncol);
    fname = [datdir dstr str1 '_50']; 
    tmp = readmatrix(fname);
    tmp = packr(tmp')';
    X_50_1(:,:,i)=tmp(n_row:n_row+1,1:ncol);
    fname = [datdir dstr str1 '_95']; 
    tmp = readmatrix(fname);
    tmp = packr(tmp')';
    X_95_1(:,:,i)=tmp(n_row:n_row+1,1:ncol);
   
    fname = [datdir dstr str5 '_5']; 
    tmp = readmatrix(fname);
    tmp = packr(tmp')';
    X_05_5(:,:,i)=tmp(n_row:n_row+1,1:ncol);
    fname = [datdir dstr str5 '_50']; 
    tmp = readmatrix(fname);
    tmp = packr(tmp')';
    X_50_5(:,:,i)=tmp(n_row:n_row+1,1:ncol);
    fname = [datdir dstr str5 '_95']; 
    tmp = readmatrix(fname);
    tmp = packr(tmp')';
    X_95_5(:,:,i)=tmp(n_row:n_row+1,1:ncol);
end
% Produce Table
outfile_name = [outdir 'um_mc_tab_' file_str '_supplementarymaterial.txt'];
fileID = fopen(outfile_name,'w');

fprintf(fileID,'Size: k = 1 \n');
fprintf(fileID,'DGP; ');
for i = 1:ncol
    if isnan(p_num(i)) == 0
      fprintf(fileID,[psrt '%4.3f'],p_num(i));
    end
    if i < ncol
        fprintf(fileID,';');
    else
        fprintf(fileID,'\n');
    end
end
for irow = 1:7
    nm = char(Design_name(irow));
    fprintf(fileID,[nm ';']);
    for icol = 1:ncol
      tmp = [X_50_1(1,icol,irow) X_05_1(1,icol,irow) X_95_1(1,icol,irow)];
      fprintf(fileID,'%4.3f (%4.3f,%4.3f)',tmp);
      if icol < ncol
        fprintf(fileID,';');
      else
        fprintf(fileID,'\n');
      end
    end
end

fprintf(fileID,'\n');
fprintf(fileID,'Size: k = 5 \n');
fprintf(fileID,'DGP; ');
for i = 1:ncol
    if isnan(p_num(i)) == 0
      fprintf(fileID,[psrt '%4.3f'],p_num(i));
    end
    if i < ncol
        fprintf(fileID,';');
    else
        fprintf(fileID,'\n');
    end
end
for irow = 1:7
    nm = char(Design_name(irow));
    fprintf(fileID,[nm ';']);
    for icol = 1:ncol
      tmp = [X_50_5(1,icol,irow) X_05_5(1,icol,irow) X_95_5(1,icol,irow)];
      fprintf(fileID,'%4.3f (%4.3f,%4.3f)',tmp);
      if icol < ncol
        fprintf(fileID,';');
      else
        fprintf(fileID,'\n');
      end
    end
end

fprintf(fileID,'\n');
fprintf(fileID,'Average Length: k = 1 \n');
fprintf(fileID,'DGP; ');
for i = 1:ncol
    if isnan(p_num(i)) == 0
      fprintf(fileID,[psrt '%4.3f'],p_num(i));
    end
    if i < ncol
        fprintf(fileID,';');
    else
        fprintf(fileID,'\n');
    end
end
for irow = 1:7
    nm = char(Design_name(irow));
    fprintf(fileID,[nm ';']);
    for icol = 1:ncol
      tmp = [X_50_1(2,icol,irow) X_05_1(2,icol,irow) X_95_1(2,icol,irow)];
      fprintf(fileID,'%4.3f (%4.3f,%4.3f)',tmp);
      if icol < ncol
        fprintf(fileID,';');
      else
        fprintf(fileID,'\n');
      end
    end
end

fprintf(fileID,'\n');
fprintf(fileID,'Average Length: k = 5 \n');
fprintf(fileID,'DGP; ');
for i = 1:ncol
    if isnan(p_num(i)) == 0
      fprintf(fileID,[psrt '%4.3f'],p_num(i));
    end
    if i < ncol
        fprintf(fileID,';');
    else
        fprintf(fileID,'\n');
    end
end
for irow = 1:7
    nm = char(Design_name(irow));
    fprintf(fileID,[nm ';']);
    for icol = 1:ncol
      tmp = [X_50_5(2,icol,irow) X_05_5(2,icol,irow) X_95_5(2,icol,irow)];
      fprintf(fileID,'%4.3f (%4.3f,%4.3f)',tmp);
      if icol < ncol
        fprintf(fileID,';');
      else
        fprintf(fileID,'\n');
      end
    end
end
      