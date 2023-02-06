clear all
input_file = '4flakes_posneg_full_style.lmp';
% input_file = '2_atoms_data.init.Shi.mol'
d = importdata(input_file,' ',9 ).data;
d(:,1)=1:length(d(:,1));
d = sortrows(d,2);
tot_mol = max(d(:,2))
tot_type = max(d(:,3))

inp = input('is the data file written in full type atom style? (y/n)')
if inp == 'n'
    error('please convert to full type atom style')
end

%%% reading box boundaries
fid=fopen(input_file); 
linenum = 4
Cx = textscan(fid,'%.6f %.6f %s',1,'delimiter','\n', 'headerlines',linenum-1);
fclose(fid);
fid=fopen(input_file); 
Cy = textscan(fid,'%.6f %.15f %s',1,'delimiter','\n', 'headerlines',linenum);
fclose(fid);
fid=fopen(input_file); 
Cz = textscan(fid,'%.6f %.6f %s',1,'delimiter','\n', 'headerlines',linenum+1);
fclose(fid);


xlo = Cx{1,1};xhi = Cx{1,2};
ylo = Cy{1,1};yhi = Cy{1,2};
zlo = Cz{1,1};zhi = Cz{1,2};

Lx = xhi-xlo; Ly = yhi-ylo; Lz = zhi-zlo; 


% chi = [5.6878 5.8678 5.8678 5.8678] %ELECTRONEGATIVITY
% J_AA   = [14 14 14 14] %% SELF COULOMB POTENTIAL

% chi = [3 1]
% J_AA = [10 10]

chi = input('Input the chi values for all the atom types: ')
J_AA   = input('Input the J_AA values for all the atom types: ')

if length(chi)~=tot_type || length(J_AA)~=tot_type
    error('not all chis or Js set')
end

d = [d zeros(length(d(:,1)),2)]

for i = 1:tot_type   
type_n(i) = sum(d(:,3)==i); %%% counting total atoms in each type note

% d(d(:,3)==i,8)=chi(i); %% placing chi and J_AA in the atom-data
% d(d(:,3)==i,9)=J_AA(i);

%%%% Matching with Prof. Shi
d(d(:,3)==i,8)=chi(i).*14.4; %% placing chi and J_AA in the atom-data
d(d(:,3)==i,9)=J_AA(i).*14.4;
end


for i = 1:tot_mol   
mol_n(i) = sum(d(:,2)==i); %%% counting atoms in each molecule. Note: d is sorted by mole-ID
tot_q_on_mol(i) = sum(d(d(:,2)==i,4))  %%  counting total charge in each molecule
end


Xrow = zeros(length(d(:,1)),length(d(:,1)));
Xcol = zeros(length(d(:,1)),length(d(:,1)));
Yrow = zeros(length(d(:,1)),length(d(:,1)));
Ycol = zeros(length(d(:,1)),length(d(:,1)));
Zrow = zeros(length(d(:,1)),length(d(:,1)));
Zcol = zeros(length(d(:,1)),length(d(:,1)));

for i = 1:length(Xrow)
Xrow(i,:) = d(:,5)';
Yrow(i,:) = d(:,6)';
Zrow(i,:) = d(:,7)';
end

for i = 1:length(Xrow)
Xcol(:,i) = d(:,5);
Ycol(:,i) = d(:,6);
Zcol(:,i) = d(:,7);
end

delX = Xrow-Xcol;delY = Yrow-Ycol;delZ = Zrow-Zcol;

%%%% Implementing boundary condition

delX(delX(:,:)>(Lx/2))=delX(delX(:,:)>(Lx/2))-Lx;
delY(delY(:,:)>(Ly/2))=delY(delY(:,:)>(Ly/2))-Ly;
delZ(delZ(:,:)>(Lz/2))=delZ(delZ(:,:)>(Lz/2))-Lz;

delX(delX(:,:)<(-Lx/2))=delX(delX(:,:)<(-Lx/2))+Lx;
delY(delY(:,:)<(-Ly/2))=delY(delY(:,:)<(-Ly/2))+Ly;
delZ(delZ(:,:)<(-Lz/2))=delZ(delZ(:,:)<(-Lz/2))+Lz;
%%% generating matrix of all distances between all atoms [r11 r12...; r21 r22 ... ;.... rnn]
del2 = delX.^2 + (delY.^2 ) +delZ.^2;
del = del2.^0.5;  %% all distances (r)

k = 14.4./del;

for i = 1:length(k)
    k(i,i) = d(i,9);
end

stt = 2;
endd = 0;

head = stt-1;
chi_full = d(:,8);
% chi_full = [3 1];
EqR = []; EqL=[];

%%%%%%%%%%%%%%%%%%%%%%%%% n-1 equations for each molecules
for j = 1:tot_mol      
    endd = endd+mol_n(j)
    EqR = [EqR;chi_full(stt:endd)-chi_full(head)];
    EqL = [EqL;k(head,:)-k(stt:endd,:)];
    stt = endd+2;  head = stt-1;
end
%%%%%%%%%%%%%%%%%%%%%%%%% total charge conservation on each molecule (1 more eq for each molecule)
QR = tot_q_on_mol';
% QL =[0;0];
% QL =  [566 0 0 -566]'; %% for testing
QL = zeros(tot_mol,length(d(:,1)));

stt = 1;
endd = 0;
for i = 1:tot_mol
    endd = endd+mol_n(i);
    QL(i,stt:endd)=1;
    stt = endd+1;
end

EqR = [EqR;QR]; EqL = [EqL;QL];

Q = EqL\EqR; %%% EqL is the C matrix, EqR is the -D matrix

% L = []; R =[];
% for j = 1:tot_mol
% L = [L;EqR{:,j}];
% R = [R;EqL{:,j}];
% end
tot_struc = [d(:,1:3) Q d(:,5:7)]; 

% fname = '4flakes_chrg_dist_matlab_100m1.lmp'
 fname= join([input_file,'_qeq.lmp'])
  delete(fname);
  
  fileID = fopen(fname,'w');
  fprintf(fileID,'# LAMMPS data file written by MATLAB\n');
  fprintf(fileID,'%d atoms\n',length(tot_struc(:,1)));
  fprintf(fileID,'%d atom types\n',max(tot_struc(:,2)));
  fprintf(fileID,'%f %f xlo xhi\n',xlo, xhi);
  fprintf(fileID,'%f %f ylo yhi\n',ylo, yhi);
  fprintf(fileID,'%f %f zlo zhi\n\n',zlo, zhi);
  fprintf(fileID,'Atoms\n');
  fclose(fileID);
  writematrix(tot_struc,fname, 'FileType','text','Delimiter','space','WriteMode','append');
  

for i = 1:tot_mol
tot_q_on_mol_after_eq(i) = sum(tot_struc(tot_struc(:,2)==i,4))
end

